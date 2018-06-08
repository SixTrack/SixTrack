#include "FlukaIOServer_private.h"
#include "Connection.h"
#include "FlukaIO.h"

/**
 * Create a server object
 * @param port Port where the server will be listening (will be randomly assigned by system if 0)
 * @return port number if success, -1 if error
 */
flukaio_server_t* flukaio_server_create()
{
	flukaio_server_t *server = (flukaio_server_t*)malloc(sizeof(flukaio_server_t));
	if (server) {
		server->port = -1;
		server->fd = -1;
		server->start = &netio_server_start;
		server->accept = &netio_server_accept;
		server->shutdown = &netio_server_shutdown;
		server->connection_create = &connection_create;
	}
	return server;
}

/**
 * Starts a server listening on port
 * @param port Port where the server will be listening (will be randomly assigned by system if 0)
 * @return port number if success, -1 if error
 */
int flukaio_server_start(flukaio_server_t* server, unsigned int port)
{
	if (server->port < 0) {
		server->port = port;
		return server->start(&server->port, &server->fd);
	} else {
		return -1;
	}
}

/**
 * Waits a new connection
 * @return Connection object with incoming connection, NULL if failure
 */
flukaio_connection_t *flukaio_server_accept(flukaio_server_t* server)
{
	int sockfd = server->accept(server->fd);
	if (sockfd < 0) return NULL;

	flukaio_connection_t *conn = server->connection_create(sockfd);

	if (flukaio_handshake_server(conn) < 0) {
		// Handshake failed
		errno = EPROTO;
		flukaio_disconnect(conn);
		return NULL;
	}

	return conn;
}

/**
 * Close and shutdown the server
 * Frees the server
 * @return 0 if success, -1 if the server was not initialized
 */
int flukaio_server_shutdown(flukaio_server_t* server)
{
	if (server) {

		if (server->fd != -1)
			server->shutdown(server->fd);

		free(server);
	}
	return 0;
}

/**
 * @return Port assigned to server
 */
int flukaio_server_get_port(flukaio_server_t* server) {
	if (server) return server->port;
	return -1;
}

