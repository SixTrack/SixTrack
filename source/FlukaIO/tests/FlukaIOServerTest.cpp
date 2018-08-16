#include "FlukaIO.h"
#include "FlukaIOServer.h"
#include "fakes/FakeNetIO.h"
#include "fakes/FakeFlukaIO.h"
#include "fakes/FakeFlukaIOServer.h"
#include "fakes/FakeFlukaIOHandshake.h"

#include "CommonTest.h"

static flukaio_server_t* server;

TEST_GROUP(FlukaIOServer)
{

	int n;

	void setup()
	{
		server = fakeflukaio_server_create();
		if (!server) {
			FAIL("Server should be valid here");
		}
	}

	void teardown()
	{
		if (server) {
			n = flukaio_server_shutdown(server);
			server = NULL;
			CHECK(n != -1);
		}
	}

};

TEST(FlukaIOServer, flukaio_server_create)
{
	LONGS_EQUAL(-1, server->port);
	LONGS_EQUAL(-1, server->fd);
}

TEST(FlukaIOServer, flukaio_server_start_random_port)
{
	int port = 234;
	fakenetio_set_next_port(port);

	n = flukaio_server_start(server, 0);

	LONGS_EQUAL(port, n);
}

TEST(FlukaIOServer, flukaio_server_start_given_port)
{
	int port = 236;

	n = flukaio_server_start(server, port);

	LONGS_EQUAL(port, n);
}

TEST(FlukaIOServer, flukaio_server_get_port_stopped)
{
	n = flukaio_server_get_port(server);

	LONGS_EQUAL(-1, n);
}

TEST(FlukaIOServer, flukaio_server_get_port)
{
	int port = 235;
	fakenetio_set_next_port(port);

	flukaio_server_start(server, 0);
	n = flukaio_server_get_port(server);

	LONGS_EQUAL(port, n);
}

TEST(FlukaIOServer, flukaio_server_start_twice)
{
	n = flukaio_server_start(server, 0);
	CHECK(n != -1);

	n = flukaio_server_start(server, 0);
	CHECK(n == -1);
}

TEST(FlukaIOServer, flukaio_server_accept)
{
	flukaio_server_start(server, 0);
	int hsk_count = fakeflukaio_handshake_server_called;

	flukaio_connection_t* conn = flukaio_server_accept(server);
	CHECK(conn != NULL);

	CHECK(fakeflukaio_handshake_server_called > hsk_count);
	flukaio_disconnect(conn);
}

TEST(FlukaIOServer, flukaio_server_shutdown_stops_the_server_if_initialized)
{
	flukaio_server_start(server, 0);
}

