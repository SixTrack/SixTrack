#include "FakeFlukaIOHandshake.h"

int (*flukaio_handshake_client)(flukaio_connection_t * conn) = &fakeflukaio_handshake_client;
int (*flukaio_handshake_server)(flukaio_connection_t * conn) = &fakeflukaio_handshake_server;
int fakeflukaio_handshake_client_called = 0;
int fakeflukaio_handshake_server_called = 0;

int fakeflukaio_handshake_client(flukaio_connection_t * conn) {
	fakeflukaio_handshake_client_called++;
	return 1;
}

int fakeflukaio_handshake_server(flukaio_connection_t * conn) {
	fakeflukaio_handshake_server_called++;
	return 1;
}
