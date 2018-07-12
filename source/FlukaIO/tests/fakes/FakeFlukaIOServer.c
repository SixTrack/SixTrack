#include "FakeFlukaIOServer.h"

flukaio_server_t* fakeflukaio_server_create() {
	flukaio_server_t *server = flukaio_server_create();
	if (server) {
		server->start = &fakenetio_server_start;
		server->accept = &fakenetio_server_accept;
		server->shutdown = &fakenetio_server_shutdown;
		server->connection_create = &fakeconnection_create;
	}
	return server;
}
