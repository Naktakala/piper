server = chi.NetworkServer.Create({
  port_number = 49468,
  --verbose = true
})

alive = true -- Variable normally set to false by HTTP client
while (alive) do
  chiNetworkServerSynchronize(server)
end

chiNetworkServerShutdown(server)