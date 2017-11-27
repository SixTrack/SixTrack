#include <stdio.h>
#include <stdint.h>
#include <stddef.h>
#include <math.h>
#include <inttypes.h>
#include <string.h>

#include "Message.h"
#include "Connection.h"
#include "FlukaIOServer.h"
#include "fakes/FakeConnection.h"
#include "fakes/FakeFortranFlukaIO.h"
#include "fakes/FakeNetIO.h"
#include "fakes/FakeFlukaIO.h"
#include "fakes/FakeFlukaIOServer.h"
#include "fakes/FakeFlukaIOHandshake.h"

#include "CommonTest.h"

TEST_GROUP(FortranFlukaIO)
{

	flukaio_connection_t *conn;
	int cid;

	particle_info_t result;
	particle_info_t expected;
	flukaio_message_t msg;

	char buf[1024];

	ssize_t n;

	void setup()
	{
		memset(&result, 0x00, sizeof(result));
		memset(&expected, 0x00, sizeof(expected));
		memset(&msg, 0x00, sizeof(msg));

		expected.id =   456;
		expected.gen = 7;
		expected.weight = 1.0;
		expected.x =   1.11e2;
		expected.y =   2.22e2;
		expected.z =   3.33e3;
		expected.tx =   11e-1;
		expected.ty =   22e-1;
		expected.tz =   33e-1;
		expected.aa =   3;
		expected.zz =   4;
		expected.m =    450e0;
		expected.pc =  4520e-1;
		expected.t =    5e-3;

		// create connection
		conn = fakeconnection_create(-1);
		ntinit_();
		cid = store_connection(conn);

	}

	void teardown()
	{
		// Check no messages are left
		n = ntrecv_(
				&cid,
				NULL,
				NULL, NULL,
				NULL,
				NULL, NULL, NULL,
				NULL, NULL, NULL,
				NULL, NULL,
				NULL, NULL,
				NULL
				);
		LONGS_EQUAL(-1, n);

		// Cleanup
		ntend_(&cid);
		conn = NULL;
	}

	void PARTICLES_EQUAL(const particle_info_t *expected, const particle_info_t *result)
	{
		LONGS_EQUAL(expected->id, result->id);
		LONGS_EQUAL(expected->gen, result->gen);
		DOUBLES_EQUAL(expected->weight, result->weight, DOUBLE_TOL);
		DOUBLES_EQUAL(expected->x, result->x, DOUBLE_TOL);
		DOUBLES_EQUAL(expected->y, result->y, DOUBLE_TOL);
		DOUBLES_EQUAL(expected->z, result->z, DOUBLE_TOL);
		DOUBLES_EQUAL(expected->tx, result->tx, DOUBLE_TOL);
		DOUBLES_EQUAL(expected->ty, result->ty, DOUBLE_TOL);
		DOUBLES_EQUAL(expected->tz, result->tz, DOUBLE_TOL);
		LONGS_EQUAL(expected->aa, result->aa);
		LONGS_EQUAL(expected->zz, result->zz);
		DOUBLES_EQUAL(expected->m, result->m, DOUBLE_TOL);
		DOUBLES_EQUAL(expected->pc, result->pc, DOUBLE_TOL);
		DOUBLES_EQUAL(expected->t, result->t, DOUBLE_TOL);
	}
	void MESSAGE_HEADER_EQUAL(flukaio_message_t *msg, flukaio_message_type_t type, size_t len) {
		LONGS_EQUAL(type, msg->type);
		LONGS_EQUAL(MSG_HEADER_LEN + len, msg->size);
	}

};

/* Async Read operations */
TEST(FortranFlukaIO, read_message_empty_buffer)
{
	uint8_t type;
	n = ntrecv_(
			&cid,
			&type,
			&result.id, &result.gen,
			&result.weight,
			&result.x, &result.y, &result.z,
			&result.tx, &result.ty, &result.tz,
			&result.aa, &result.zz,
			&result.m, &result.pc, &result.t);

	LONGS_EQUAL(-1, n);
}

TEST(FortranFlukaIO, read_NULL_message_particle)
{
	// Read in NULL buffer skips next message
	fakeflukaio_insert_particle_message(conn, &expected);

	n = ntrecv_(
			&cid,
			NULL,
			NULL, NULL,
			NULL,
			NULL, NULL, NULL,
			NULL, NULL, NULL,
			NULL, NULL,
			NULL, NULL,
			NULL
			);
	LONGS_EQUAL(MSG_PART_LEN, n);
}


TEST(FortranFlukaIO, read_message_particle)
{
	fakeflukaio_insert_particle_message(conn, &expected);

	uint8_t type;
	n = ntrecv_(
			&cid,
			&type,
			&result.id, &result.gen,
			&result.weight,
			&result.x, &result.y, &result.z,
			&result.tx, &result.ty, &result.tz,
			&result.aa, &result.zz,
			&result.m, &result.pc, &result.t);

	CHECK(n != -1);
	LONGS_EQUAL(N_PART, type);
	PARTICLES_EQUAL(&expected, &result);
}

TEST(FortranFlukaIO, read_message_eob)
{
	fakeflukaio_insert_eob_message(conn);

	uint8_t type;
	n = ntrecv_(
			&cid,
			&type,
			&result.id, &result.gen,
			&result.weight,
			&result.x, &result.y, &result.z,
			&result.tx, &result.ty, &result.tz,
			&result.aa, &result.zz,
			&result.m, &result.pc, &result.t);
	CHECK(n != -1);

	LONGS_EQUAL(N_EOB, type);
}

TEST(FortranFlukaIO, read_message_insertion_point)
{
	uint32_t turn = 5694;
	uint16_t ipt = 8;
	fakeflukaio_insert_ipt_message(conn, turn, ipt);

	uint8_t type;
	n = ntrecv_(
			&cid,
			&type,
			&result.id, &result.gen,
			&result.weight,
			&result.x, &result.y, &result.z,
			&result.tx, &result.ty, &result.tz,
			&result.aa, &result.zz,
			&result.m, &result.pc, &result.t);
	CHECK(n != -1);

	LONGS_EQUAL(N_IPT, type);
	LONGS_EQUAL(turn, result.gen);
	LONGS_EQUAL(ipt, result.id);
}

TEST(FortranFlukaIO, ntrecv_with_eoc)
{
	fakeflukaio_insert_eoc_message(conn);

	uint8_t type;
	n = ntrecv_(
			&cid,
			&type,
			&result.id, &result.gen,
			&result.weight,
			&result.x, &result.y, &result.z,
			&result.tx, &result.ty, &result.tz,
			&result.aa, &result.zz,
			&result.m, &result.pc, &result.t);
	CHECK(n != -1);

	LONGS_EQUAL(N_END, type);
}

TEST(FortranFlukaIO, ntrecv_with_particle_and_end_of_turn)
{
	fakeflukaio_insert_particle_message(conn, &expected);
	fakeflukaio_insert_eob_message(conn);

	uint8_t type;
	n = ntrecv_(
			&cid,
			&type,
			&result.id, &result.gen,
			&result.weight,
			&result.x, &result.y, &result.z,
			&result.tx, &result.ty, &result.tz,
			&result.aa, &result.zz,
			&result.m, &result.pc, &result.t);
	CHECK(n != -1);
	LONGS_EQUAL(N_PART, type);
	PARTICLES_EQUAL(&expected, &result);

	n = ntrecv_(
			&cid,
			&type,
			&result.id, &result.gen,
			&result.weight,
			&result.x, &result.y, &result.z,
			&result.tx, &result.ty, &result.tz,
			&result.aa, &result.zz,
			&result.m, &result.pc, &result.t);

	CHECK(n != -1);
	LONGS_EQUAL(N_EOB, type);
}

TEST(FortranFlukaIO, ntrecv_with_other_message)
{
	fakeflukaio_insert_unknown_message(conn);
	uint8_t type;
	n = ntrecv_(
			&cid,
			&type,
			&result.id, &result.gen,
			&result.weight,
			&result.x, &result.y, &result.z,
			&result.tx, &result.ty, &result.tz,
			&result.aa, &result.zz,
			&result.m, &result.pc, &result.t);

	LONGS_EQUAL(-1, n);
}

/* Send operations */
// Send end of batch
TEST(FortranFlukaIO, ntsend_eob)
{
	n = ntsendeob_(&cid);
	LONGS_EQUAL(MSG_EOB_LEN, fakenetio_get_last_sent_size());

	fakeflukaio_get_last_sent_message(conn, &msg);
	MESSAGE_HEADER_EQUAL(&msg, N_EOB, 0);

}

// Send iot
TEST(FortranFlukaIO, ntsend_insertion_point)
{
	uint32_t turn = 1354;
	uint16_t ipt = 7;
	n = ntsendipt_(&cid, &turn, &ipt);

	connection_flush(conn);
	LONGS_EQUAL(MSG_IPT_LEN, fakenetio_get_last_sent_size());

	fakeflukaio_get_last_sent_message(conn, &msg);
	MESSAGE_HEADER_EQUAL(&msg, N_IPT, sizeof(flukaio_ipt_data_t));
	LONGS_EQUAL(turn, msg.data.ipt.turn);
	LONGS_EQUAL(ipt, msg.data.ipt.ipt);

}

// Send Particle
TEST(FortranFlukaIO, send_particle)
{
	n = ntsendp_(
			&cid,
			&expected.id, &expected.gen,
			&expected.weight,
			&expected.x, &expected.y, &expected.z,
			&expected.tx, &expected.ty, &expected.tz,
			&expected.aa, &expected.zz,
			&expected.m, &expected.pc,
			&expected.t
			);

	connection_flush(conn);
	LONGS_EQUAL(MSG_PART_LEN, fakenetio_get_last_sent_size());

	fakeflukaio_get_last_sent_message(conn, &msg);
	MESSAGE_HEADER_EQUAL(&msg, N_PART, sizeof(expected));
	PARTICLES_EQUAL(&expected, &msg.data.particle);
}

// Send Particle n times
TEST(FortranFlukaIO, send_particle_n_times)
{
	int times = 25;
	for(int i = 0; i < times; ++i) {
		expected.id = i;
		n = ntsendp_(
				&cid,
				&expected.id, &expected.gen,
				&expected.weight,
				&expected.x, &expected.y, &expected.z,
				&expected.tx, &expected.ty, &expected.tz,
				&expected.aa, &expected.zz,
				&expected.m, &expected.pc,
				&expected.t
				);
		LONGS_EQUAL(MSG_PART_LEN, n);
	}

	connection_flush(conn);
	LONGS_EQUAL(MSG_PART_LEN * times, fakenetio_get_last_sent_size());

	for(int i = 0; i < times; ++i) {
		expected.id = i;
		fakeflukaio_get_last_sent_message(conn, &msg);
		MESSAGE_HEADER_EQUAL(&msg, N_PART, sizeof(expected));
		PARTICLES_EQUAL(&expected, &msg.data.particle);
	}
}


// End Computation
TEST(FortranFlukaIO, send_end_comp)
{
	n = ntsendeoc_(&cid);
	LONGS_EQUAL(MSG_EOC_LEN, fakenetio_get_last_sent_size());

	fakeflukaio_get_last_sent_message(conn, &msg);
	MESSAGE_HEADER_EQUAL(&msg, N_END, 0);
}

/* Sync Read operations */

TEST(FortranFlukaIO, ntwait_with_empty_buffer)
{
	uint8_t type;
	n = ntwait_(
			&cid,
			&type,
			&result.id, &result.gen,
			&result.weight,
			&result.x, &result.y, &result.z,
			&result.tx, &result.ty, &result.tz,
			&result.aa, &result.zz,
			&result.m, &result.pc,
			&result.t);

	LONGS_EQUAL(-1, n);
}

TEST(FortranFlukaIO, ntwait_with_one_particle)
{
	fakeflukaio_insert_particle_message(conn, &expected);

	uint8_t type;
	n = ntwait_(
			&cid,
			&type,
			&result.id, &result.gen,
			&result.weight,
			&result.x, &result.y, &result.z,
			&result.tx, &result.ty, &result.tz,
			&result.aa, &result.zz,
			&result.m, &result.pc,
			&result.t);

	LONGS_EQUAL(MSG_PART_LEN, n);
	LONGS_EQUAL(N_PART, type);
	PARTICLES_EQUAL(&expected, &result);
}

TEST(FortranFlukaIO, ntwait_with_eob)
{
	fakeflukaio_insert_eob_message(conn);

	uint8_t type;
	n = ntwait_(
			&cid,
			&type,
			&result.id, &result.gen,
			&result.weight,
			&result.x, &result.y, &result.z,
			&result.tx, &result.ty, &result.tz,
			&result.aa, &result.zz,
			&result.m, &result.pc,
			&result.t);

	LONGS_EQUAL(MSG_EOB_LEN, n);
	LONGS_EQUAL(N_EOB, type);
}

TEST(FortranFlukaIO, ntwait_with_ipt)
{
	fakeflukaio_insert_ipt_message(conn, 12, 9);

	uint8_t type;
	n = ntwait_(
			&cid,
			&type,
			&result.id, &result.gen,
			&result.weight,
			&result.x, &result.y, &result.z,
			&result.tx, &result.ty, &result.tz,
			&result.aa, &result.zz,
			&result.m, &result.pc,
			&result.t);

	LONGS_EQUAL(MSG_IPT_LEN, n);
	LONGS_EQUAL(N_IPT, type);
	LONGS_EQUAL(12, result.gen);
	LONGS_EQUAL(9, result.id);
}

TEST(FortranFlukaIO, ntwait_with_eoc)
{
	fakeflukaio_insert_eoc_message(conn);

	uint8_t type;
	n = ntwait_(
			&cid,
			&type,
			&result.id, &result.gen,
			&result.weight,
			&result.x, &result.y, &result.z,
			&result.tx, &result.ty, &result.tz,
			&result.aa, &result.zz,
			&result.m, &result.pc,
			&result.t);

	LONGS_EQUAL(MSG_EOC_LEN, n);
	LONGS_EQUAL(N_END, type);
	LONGS_EQUAL(N_END, type);
}

TEST(FortranFlukaIO, ntwait_with_ipt_part_and_end_of_turn)
{
	fakeflukaio_insert_ipt_message(conn, 66, 3);
	fakeflukaio_insert_particle_message(conn, &expected);
	fakeflukaio_insert_eob_message(conn);

	uint8_t type;

	n = ntwait_(
			&cid,
			&type,
			&result.id, &result.gen,
			&result.weight,
			&result.x, &result.y, &result.z,
			&result.tx, &result.ty, &result.tz,
			&result.aa, &result.zz,
			&result.m, &result.pc,
			&result.t);
	LONGS_EQUAL(MSG_IPT_LEN, n);
	LONGS_EQUAL(N_IPT, type);
	LONGS_EQUAL(3, result.id);
	LONGS_EQUAL(66, result.gen);

	n = ntwait_(
			&cid,
			&type,
			&result.id, &result.gen,
			&result.weight,
			&result.x, &result.y, &result.z,
			&result.tx, &result.ty, &result.tz,
			&result.aa, &result.zz,
			&result.m, &result.pc,
			&result.t);
	LONGS_EQUAL(MSG_PART_LEN, n);
	LONGS_EQUAL(N_PART, type);
	PARTICLES_EQUAL(&expected, &result);

	n = ntwait_(
			&cid,
			&type,
			&result.id, &result.gen,
			&result.weight,
			&result.x, &result.y, &result.z,
			&result.tx, &result.ty, &result.tz,
			&result.aa, &result.zz,
			&result.m, &result.pc,
			&result.t);

	LONGS_EQUAL(MSG_EOB_LEN, n);
	LONGS_EQUAL(N_EOB, type);
}

TEST(FortranFlukaIO, ntwait_with_other_message)
{
	fakeflukaio_insert_unknown_message(conn);
	uint8_t type;
	n = ntwait_(
			&cid,
			&type,
			&result.id, &result.gen,
			&result.weight,
			&result.x, &result.y, &result.z,
			&result.tx, &result.ty, &result.tz,
			&result.aa, &result.zz,
			&result.m, &result.pc,
			&result.t);

	LONGS_EQUAL(-1, n);
}

static flukaio_server_t* server;

TEST_GROUP(FortranFlukaIOServer)
{
	ssize_t n;
	int serverid;


	void setup()
	{
		server = fakeflukaio_server_create();
		if (!server) {
			FAIL("Server should be valid here");
		}
		serverid = store_server(server);
		if (serverid < 0) {
			FAIL("Serverid should be valid here");
		}
	}

	void teardown()
	{
		if (serverid > -1) {
			ntshdwn_(&serverid);
		}
	}
};
/* FlukaIOServer interface */
TEST(FortranFlukaIOServer, flukaio_server_create)
{
	LONGS_EQUAL(-1, ntgetport_(&serverid));
}

TEST(FortranFlukaIOServer, flukaio_server_start_random_port)
{
	int port = 234;
	int sport = 0;
	fakenetio_set_next_port(port);

	n = ntstart_(&serverid, &sport);

	LONGS_EQUAL(port, n);
}

TEST(FortranFlukaIOServer, flukaio_server_start_given_port)
{
	int port = 236;

	n = ntstart_(&serverid, &port);

	LONGS_EQUAL(port, n);
}

TEST(FortranFlukaIOServer, flukaio_server_get_port_stopped)
{
	n = ntgetport_(&serverid);

	LONGS_EQUAL(-1, n);
}

TEST(FortranFlukaIOServer, flukaio_server_get_port)
{
	int port = 235;
	int sport = 0;
	fakenetio_set_next_port(port);

	ntstart_(&serverid, &sport);
	n = ntgetport_(&serverid);

	LONGS_EQUAL(port, n);
}

TEST(FortranFlukaIOServer, flukaio_server_start_twice)
{
	int sport = 0;
	n = ntstart_(&serverid, &sport);
	CHECK(n != -1);

	sport = 0;
	n = ntstart_(&serverid, &sport);
	CHECK(n == -1);
}

TEST(FortranFlukaIOServer, flukaio_server_accept)
{
	int sport = 0;
	int cid;
	n = ntstart_(&serverid, &sport);
	int hsk_count = fakeflukaio_handshake_server_called;

	cid = ntaccept_(&serverid);
	CHECK(cid != -1);

	CHECK(fakeflukaio_handshake_server_called > hsk_count);
	ntend_(&cid);
}

TEST(FortranFlukaIOServer, flukaio_server_shutdown_fails_if_not_initialized)
{
	n = ntshdwn_(&serverid);
	serverid = -1;

	CHECK(n != -1);
}

TEST(FortranFlukaIOServer, flukaio_server_shutdown_stops_the_server_if_initialized)
{
	int sport = 0;
	n = ntstart_(&serverid, &sport);
	CHECK(n != -1);

	n = ntshdwn_(&serverid);
	serverid = -1;

	CHECK(n != -1);
}

TEST(FortranFlukaIOServer, fortran_string_conversion_spaces_at_end)
{
	const char *fortran_string = "localhost                        ";
	char *new_string = NULL;
	new_string = create_fortran_string(fortran_string, strlen(fortran_string));

	CHECK(new_string != NULL);
	LONGS_EQUAL(9, strlen(new_string));
	CHECK(strcmp(new_string, "localhost") == 0);
	free(new_string);
	new_string = NULL;
}

TEST(FortranFlukaIOServer, fortran_string_conversion_empty_string)
{
	const char *fortran_string = "";
	char *new_string = NULL;
	new_string = create_fortran_string(fortran_string, strlen(fortran_string));

	CHECK(new_string != NULL);
	LONGS_EQUAL(0, strlen(new_string));
	CHECK(strcmp(new_string, "") == 0);
	free(new_string);
	new_string = NULL;
}

