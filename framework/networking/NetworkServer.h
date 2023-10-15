#ifndef CHITECH_NETWORKSERVER_H
#define CHITECH_NETWORKSERVER_H

#include "ChiObject.h"

#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>

#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>

namespace chi
{

class NetworkServer : public ChiObject
{
public:
  static InputParameters GetInputParameters();
  explicit NetworkServer(const InputParameters& params);

  /**Lua wrapped method called in a loop such that a network message can be
  * received by location zero, broadcast to other processors, the responses
  * collected and send back via to the client.*/
  void Synchronize();

  /**Lua-wrapped method to shutdown the listening server.*/
  void ShutdownServer();

  /**Destructor that will call ShutdownServer.*/
  ~NetworkServer();

protected:
  /**Meant to be run in a thread. This method essentially implements a
  * listening server. It is fairly cheap since it waits for connections in a
  * blocking fashion.*/
  void Listen();

  /**Assuming the buffered message is standard HTML format,
   * https://www.tutorialspoint.com/http/http_requests.htm, this method
   * firstly looks at the first line of the message. If a `POST` method used
   * then it will split the message into lines and collect the lines of data
   * passed the header (first blank line).*/
  std::vector<std::string>
  ProcessMessage(const std::string& buffered_message);

  /**This method executes a console command and cannot allow return values.
  * This allows lua commands to be assigned better.*/
  static ParameterBlock ConsoleDoString(const std::string& the_string);

  /**Builds a standard response when no data supplied in the request.*/
  static std::string BuildStandardResponse();

  // Input parameters
  const bool verbose_;
  const int port_number_;
  const std::string ip_address_;
  const int connection_que_size_;

  // Runtime parameters
  std::string incoming_mode_;
  std::vector<std::string> incoming_message_lines_;
  std::string outgoing_message_;

  int server_socket_ = -1;
  int client_socket_ = -1;

  bool quit_thread_flag_ = false;
  std::thread listening_thread_;
  std::mutex incoming_message_queue_mutex_;
  std::mutex outgoing_message_queue_mutex_;

  std::condition_variable thread_condition_variable_;
  bool outgoing_message_ready_ = false;
};

} // namespace chi

#endif // CHITECH_NETWORKSERVER_H
