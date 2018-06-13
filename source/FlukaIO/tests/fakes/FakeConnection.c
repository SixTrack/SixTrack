#include "FakeConnection.h"
#include "FakeNetIO.h"
#include "FakeFlukaIO.h"

flukaio_connection_t* fakeconnection_create(int fd) {
	flukaio_connection_t *conn = connection_create(fd);
	if (!conn) return NULL;

	// Fake connection
	conn->connect = &fakenetio_connect;
	conn->read = &fakenetio_read;
	conn->write = &fakenetio_write;
	conn->can_read = &fakenetio_can_read;
	conn->can_write = &fakenetio_can_write;
	conn->set_nonblocking = &fakenetio_set_nonblocking;
	conn->set_nodelay = &fakenetio_set_nodelay;

	return conn;
}
