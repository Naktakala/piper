#include "NetworkServer.h"

#include "ChiObjectFactory.h"
#include "utils/chi_utils.h"
#include "chi_log.h"
#include "console/chi_console.h"
#include "chi_lua.h"
#include "utils/chi_timer.h"

#include <unistd.h>
#include <functional>
#include <fcntl.h>

#define CHITECH_NETWORKSERVER_RECV_BUFFER_SIZE 4096

namespace chi
{

RegisterChiObject(chi, NetworkServer);

chi::InputParameters NetworkServer::GetInputParameters()
{
  InputParameters params = ChiObject::GetInputParameters();

  params.SetGeneralDescription("Networking server");

  params.AddRequiredParameter<int>("port_number", "Port number");
  params.AddOptionalParameter(
    "ip_address", "0.0.0.0", "IP address for the server");

  params.AddOptionalParameter(
    "connection_que_size", 5, "Queue size for incoming connections");

  params.AddOptionalParameter("verbose",
                              false,
                              "Controls verbosity of the server. If true will"
                              "display messages received.");

  return params;
}

/**Constructor.*/
NetworkServer::NetworkServer(const InputParameters& params)
  : ChiObject(params),
    verbose_(params.GetParamValue<bool>("verbose")),
    port_number_(params.GetParamValue<int>("port_number")),
    ip_address_(params.GetParamValue<std::string>("ip_address")),
    connection_que_size_(params.GetParamValue<int>("connection_que_size"))
{
  if (Chi::mpi.location_id == 0)
  {
    Chi::log.LogAll() << "Creating NetworkServer at " << ip_address_ << ":"
                      << port_number_;

    server_socket_ = socket(/*domain*/ AF_INET,
                            /*type*/ SOCK_STREAM,
                            /*protocol*/ 0);

    ChiLogicalErrorIf(server_socket_ < 0, "Error opening socket.");

    sockaddr_in socket_address{0, 0, 0, {0}, {'\0'}};

    socket_address.sin_family = AF_INET;
    socket_address.sin_addr.s_addr = inet_addr(ip_address_.c_str());
    socket_address.sin_port = htons(port_number_);

    int reusePort = 1;
    setsockopt(
      server_socket_, SOL_SOCKET, SO_REUSEPORT, &reusePort, sizeof(reusePort));

    ChiLogicalErrorIf(bind(server_socket_,
                           reinterpret_cast<sockaddr*>(&socket_address),
                           sizeof(socket_address)) < 0,
                      "Error binding socket");

    Chi::log.LogAll() << "Server started at " << ip_address_ << ":"
                      << port_number_;
    listening_thread_ = std::thread(
      std::bind(&chi::NetworkServer::Listen, this));
  }
}

/**Lua wrapped method called in a loop such that a network message can be
 * received by location zero, broadcast to other processors, the responses
 * collected and send back via to the client.*/
void NetworkServer::Synchronize()
{
  std::string consolidated_messages;
  incoming_message_queue_mutex_.lock();
  {
    for (const auto& message : incoming_message_lines_)
      consolidated_messages.append(message);
    incoming_message_lines_.clear();
  }
  incoming_message_queue_mutex_.unlock();

  size_t message_size = consolidated_messages.size();

  MPI_Bcast(&message_size,  // send/recv buffer
            1,              // count
            MPI_UINT64_T,   // datatype
            0,              // root
            Chi::mpi.comm); // communicator

  if (message_size > 0)
  {
    if (Chi::mpi.location_id != 0)
      consolidated_messages.assign(message_size, ' ');

    MPI_Bcast(consolidated_messages.data(),   // send/recv buffer
              static_cast<int>(message_size), // count
              MPI_CHAR,                       // datatype
              0,                              // root
              Chi::mpi.comm);                 // communicator

    auto commands = chi::StringSplit(consolidated_messages, ";");

    std::string consolidated_results;
    for (const auto& command : commands)
    {
      auto result = ConsoleDoString(command);
      result.RecursiveDumpToJSON(consolidated_results);
      consolidated_results += ";";
    }
    std::stringstream outstr;
    outstr << "HTTP/1.1 200 OK\n";
    outstr << "Content-Length: " << consolidated_results.size() << "\n\n";
    outstr << consolidated_results;

    outgoing_message_ = outstr.str();
    {
      std::lock_guard lk(outgoing_message_queue_mutex_);
      outgoing_message_ready_ = true;
    }
    thread_condition_variable_.notify_one();
  }
}

/**Lua-wrapped method to shutdown the listening server.*/
void NetworkServer::ShutdownServer()
{
  close(server_socket_);
  close(client_socket_);

  try {
    if (Chi::mpi.location_id == 0)
    {
      quit_thread_flag_ = true;
      listening_thread_.join();
    }
  }
  catch (...)
  {
    Chi::log.LogAll() << "Network thread cleanup failed";
  }

  Chi::log.LogAll() << "Server at " << ip_address_ << ":" << port_number_
                    << " stopped.";
}

/**Destructor that will call ShutdownServer.*/
NetworkServer::~NetworkServer() { ShutdownServer(); }

/**Meant to be run in a thread. This method essentially implements a
 * listening server. It is fairly cheap since it waits for connections in a
 * blocking fashion.*/
void NetworkServer::Listen()
{
  while (not quit_thread_flag_)
  {
    Chi::log.LogAll() << "Listening ";
    const int err = listen(server_socket_, connection_que_size_);
    if (err)
    {
      Chi::log.LogAllError() << "Socket listen failed";
      break ;
    }

    sockaddr_in client_socket_address{0, 0, 0, {0}, {'\0'}};
    unsigned int client_socket_len = sizeof(client_socket_address);

    client_socket_ = accept(
      server_socket_, (sockaddr*)&client_socket_address, &client_socket_len);

    const std::string client_ip_address =
      inet_ntoa(client_socket_address.sin_addr);
    const int client_port = ntohs(client_socket_address.sin_port);

    ChiLogicalErrorIf(client_socket_ < 0,
                      "Incoming connection failed from " + client_ip_address +
                        " port " + std::to_string(client_port));

    Chi::log.LogAll() << "Server connection accepted from Client " +
                           client_ip_address + ":" +
                           std::to_string(client_port);

    std::string incoming_message;
    char buffer[CHITECH_NETWORKSERVER_RECV_BUFFER_SIZE] = {0};
    long recv_size;
    bool client_connected = true;
    while (client_connected)
    {
      do
      {
        recv_size =
          read(client_socket_, buffer, CHITECH_NETWORKSERVER_RECV_BUFFER_SIZE);

        incoming_message.clear();
        for (long c = 0; c < recv_size; ++c)
          incoming_message += buffer[c];

        if (recv_size <= CHITECH_NETWORKSERVER_RECV_BUFFER_SIZE) break;
      } while (recv_size > 0);

      if (recv_size <= 0) break; // client disconnected

      incoming_message_queue_mutex_.lock();
      {
        incoming_message_lines_ = ProcessMessage(incoming_message);
      }
      incoming_message_queue_mutex_.unlock();

      if (verbose_)
        Chi::log.LogAll() << "Incoming message num lines "
                          << incoming_message_lines_.size();

      if (incoming_message_lines_.empty())
      {
        const std::string std_response = BuildStandardResponse();
        const size_t bytessent =
          write(client_socket_, std_response.c_str(), std_response.size());

        if (bytessent != std_response.size())
        {
          Chi::log.LogAllError() << "Error sending default response";
          client_connected = false;
        }
      } // if incoming_message_lines_.empty()
      else
      {
        std::unique_lock lk(outgoing_message_queue_mutex_, std::defer_lock);
        thread_condition_variable_.wait(
          lk, [&] { return outgoing_message_ready_; });

        const size_t bytessent = write(
          client_socket_, outgoing_message_.c_str(), outgoing_message_.size());

        if (bytessent != outgoing_message_.size())
        {
          Chi::log.LogAllError() << "Error sending message";
          client_connected = false;
        }

        outgoing_message_.clear();
        outgoing_message_ready_ = false;
      } // if NOT incoming_message_lines_.empty()
    }   // while client connected
    Chi::log.LogAll() << "Client " + client_ip_address + ":" +
                           std::to_string(client_port)
                      << " Disconnected";
  } // while not quit posted
}

/**Assuming the buffered message is standard HTML format,
 * https://www.tutorialspoint.com/http/http_requests.htm, this method
 * firstly looks at the first line of the message. If a `POST` method used then
 * it will split the message into lines and collect the lines of data passed
 * the header (first blank line).*/
std::vector<std::string>
NetworkServer::ProcessMessage(const std::string& buffered_message)
{
  constexpr size_t NPOS = std::string::npos;
  //============================================= Lambda definition
  /**This lambda split a string based on a delimiter. It is specially modified
   * to return empty pieces even successive delimeters are encountered.*/
  auto ModifiedStringSplit =
    [](const std::string& input, const std::string& delim /*=" "*/)
  {
    constexpr size_t NPOS = std::string::npos;
    std::vector<std::string> output;

    std::string remainder = input;
    size_t first_scope = remainder.find_first_of(delim);

    while (first_scope != NPOS)
    {
      if (first_scope != 0) output.push_back(remainder.substr(0, first_scope));

      if (remainder.substr(0, delim.size()) == delim) output.push_back("");

      remainder = remainder.substr(first_scope + delim.size(), NPOS);
      first_scope = remainder.find_first_of(delim);
    }
    output.push_back(remainder);

    return output;
  };

  //============================================= Split off first line
  // We want the standard HTML syntax,
  // see https://www.tutorialspoint.com/http/http_requests.htm
  const size_t end_line1 = buffered_message.find_first_of('\n');

  if (end_line1 == NPOS) return {};

  const std::string line1 = buffered_message.substr(0, end_line1);
  const auto words = chi::StringSplit(line1);

  if (words.size() < 2)
  {
    Chi::log.Log0Error() << "HTTP request: \"" << line1
                         << "\" cannot be processed.";
    return {};
  }

  //============================================= Determine the command
  const std::string& request_method = words[0];

  if (request_method == "POST" or request_method == "GET")
  {
    incoming_mode_ = words[0];
    //====================================== Split the message into lines
    // Standard HTML has carriage return \r AND line-feed \n
    const auto lines = ModifiedStringSplit(buffered_message, "\r\n");

    //====================================== Skip till empty line then collect
    std::vector<std::string> data_lines;
    bool data_lines_active = false;
    for (const auto& line : lines)
    {
      if (line.empty())
      {
        data_lines_active = true;
        continue;
      }

      if (data_lines_active) data_lines.push_back(line);
    }

    return data_lines;
  }
  else
    incoming_mode_ = "";
  return {};
}

/**This method executes a console command and returns "abstractly" what is at
 * the top of the stack.*/
ParameterBlock NetworkServer::ConsoleDoString(const std::string& the_string)
{
  lua_State* L = Chi::console.GetConsoleState();

  lua_settop(L, 0);
  bool error = luaL_dostring(L, the_string.c_str());
  if (error)
  {
    Chi::log.LogAll() << lua_tostring(L, -1);
    lua_pop(L, 1);
  }
  const int num_retvals = lua_gettop(L);
  if (num_retvals == 0) return ParameterBlock{};
  else
    return chi_lua::StackItemToParameterBlock(L, -1);
}

/**Builds a standard response when no data supplied in the request.*/
std::string NetworkServer::BuildStandardResponse()
{
  const std::string html_page =
    "<!DOCTYPE html><html lang=\"en\"><body><h1> HOME </h1><p> "
    "ChiTech:InvalidRequest :) </p></body></html>";
  std::ostringstream ss;
  ss << "HTTP/1.1 200 OK\nContent-Type: text/html\nContent-Length: "
     << html_page.size() << "\n\n"
     << html_page << "\n";

  return ss.str();
}

} // namespace chi