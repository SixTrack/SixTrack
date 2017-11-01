#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <stddef.h>
#include <inttypes.h>
#include <string.h>

#include "Message.h"
#include "Connection.h"
#include "fakes/FakeNetIO.h"
#include "fakes/FakeFlukaIO.h"
#include "fakes/FakeConnection.h"
#include "fakes/FakeFlukaIOHandshake.h"

#include "CommonTest.h"

TEST_GROUP(FlukaIO)
{

    flukaio_connection_t *conn;

    particle_info_t result;
    particle_info_t expected;
    flukaio_message_t msg;

    char buf[1024];

    ssize_t n;

    void setup()
    {

        // create connection
        conn = fakeconnection_create(-1);

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

    }

    void teardown()
    {
        // Check no messages are left
        n = flukaio_receive_message(conn, NULL);
        LONGS_EQUAL(-1, n);

        // Cleanup
        flukaio_disconnect(conn);
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

/* connect */
TEST(FlukaIO, connect)
{
    conn = flukaio_connect(conn, "htest", 0);
	CHECK(conn != NULL);
    CHECK(fakeflukaio_handshake_client_called == 1);
}

/* Read operations */
TEST(FlukaIO, read_message_empty_buffer)
{
    n = flukaio_receive_message(conn, &msg);

    LONGS_EQUAL(-1, n);
}

TEST(FlukaIO, read_NULL_message_particle)
{
    // Read in NULL buffer skips next message
    fakeflukaio_insert_particle_message(conn, &expected);

    n = flukaio_receive_message(conn, NULL);
    LONGS_EQUAL(MSG_PART_LEN, n);
}

TEST(FlukaIO, read_message_particle)
{
    // Read particle message
    fakeflukaio_insert_particle_message(conn, &expected);

    n = flukaio_receive_message(conn, &msg);
    CHECK(n != -1);

    MESSAGE_HEADER_EQUAL(&msg, N_PART, sizeof(expected));
    PARTICLES_EQUAL(&expected, &msg.data.particle);
}

TEST(FlukaIO, read_message_particle_twice)
{
    // Insert two messages and check they are read correctly
    fakeflukaio_insert_particle_message(conn, &expected);
    fakeflukaio_insert_particle_message(conn, &expected);

    n = flukaio_receive_message(conn, &msg);
    CHECK(n != -1);

    MESSAGE_HEADER_EQUAL(&msg, N_PART, sizeof(expected));
    PARTICLES_EQUAL(&expected, &msg.data.particle);

    memset(&msg, 0xff, sizeof(msg));
    n = flukaio_receive_message(conn, &msg);
    CHECK(n != -1);

    MESSAGE_HEADER_EQUAL(&msg, N_PART, sizeof(expected));
    PARTICLES_EQUAL(&expected, &msg.data.particle);
}

TEST(FlukaIO, read_message_eob)
{
    // Insert end of turn and check
    fakeflukaio_insert_eob_message(conn);

    n = flukaio_receive_message(conn, &msg);
    CHECK(n != -1);

    MESSAGE_HEADER_EQUAL(&msg, N_EOB, 0);
}

TEST(FlukaIO, read_message_insertion_point)
{
    // Insert end of turn and check
    uint32_t turn = 7;
    uint16_t ipt = 2;
    fakeflukaio_insert_ipt_message(conn, turn, ipt);

    n = flukaio_receive_message(conn, &msg);
    CHECK(n != -1);

    MESSAGE_HEADER_EQUAL(&msg, N_IPT, sizeof(flukaio_ipt_data_t));
    LONGS_EQUAL(ipt, msg.data.ipt.ipt);
}

TEST(FlukaIO, read_message_insertion_point_last)
{
    // Insert end of turn and check
    uint32_t turn = -1;
    uint16_t ipt = 2;
    fakeflukaio_insert_ipt_message(conn, turn, ipt);

    n = flukaio_receive_message(conn, &msg);
    CHECK(n != -1);

    MESSAGE_HEADER_EQUAL(&msg, N_IPT, sizeof(flukaio_ipt_data_t));
    LONGS_EQUAL(turn, msg.data.ipt.turn);
    LONGS_EQUAL(ipt, msg.data.ipt.ipt);
}

TEST(FlukaIO, read_message_particle_and_end_of_turn)
{
    // Insert particle and end of turn, check
    fakeflukaio_insert_particle_message(conn, &expected);
    fakeflukaio_insert_eob_message(conn);

    // read the particle
    n = flukaio_receive_message(conn, &msg);
    CHECK(n != -1);

    MESSAGE_HEADER_EQUAL(&msg, N_PART, sizeof(expected));
    PARTICLES_EQUAL(&expected, &msg.data.particle);

    // read the turn
    n = flukaio_receive_message(conn, &msg);
    CHECK(n != -1);

    MESSAGE_HEADER_EQUAL(&msg, N_EOB, 0);
}

TEST(FlukaIO, read_message_end_of_computation)
{
    // Insert end of computation and check
    fakeflukaio_insert_eoc_message(conn);

    n = flukaio_receive_message(conn, &msg);
    CHECK(n != -1);

    MESSAGE_HEADER_EQUAL(&msg, N_END, 0);
}

TEST(FlukaIO, read_message_configuration)
{
    // Insert configuration message and check
    fakeflukaio_insert_config_message(conn);

    n = flukaio_receive_message(conn, &msg);
    CHECK(n != -1);

    MESSAGE_HEADER_EQUAL(&msg, N_CONF, 0);
}

TEST(FlukaIO, read_message_other)
{
    // Insert other message and check error is read
    fakeflukaio_insert_unknown_message(conn);

    n = flukaio_receive_message(conn, &msg);
	LONGS_EQUAL(-1, n);
    CHECK(n == -1);

    MESSAGE_HEADER_EQUAL(&msg, N_ERR, 0);
}

TEST(FlukaIO, read_message_handshake)
{
    // Insert end of turn and check
    uint16_t major = 3;
    uint16_t minor = 1;
    uint32_t key   = 123;
    fakeflukaio_insert_hsk_message(conn, major, minor, key);

    n = flukaio_receive_message(conn, &msg);
    CHECK(n != -1);

    MESSAGE_HEADER_EQUAL(&msg, N_HSK, sizeof(flukaio_hsk_data_t));
    LONGS_EQUAL(major, msg.data.hsk.major);
    LONGS_EQUAL(minor, msg.data.hsk.minor);
    LONGS_EQUAL(key,   msg.data.hsk.key);
}

/* Send operations */
/* Send end of turn */
TEST(FlukaIO, send_eob)
{
    n = flukaio_send_eob(conn);
    LONGS_EQUAL(MSG_EOB_LEN, fakenetio_get_last_sent_size());

    fakeflukaio_get_last_sent_message(conn, &msg);
    MESSAGE_HEADER_EQUAL(&msg, N_EOB, 0);
}

/* Send end of insertion point */
TEST(FlukaIO, send_insertion_point)
{
    uint32_t turn = 32;
    uint16_t ipt = 55;
    n = flukaio_send_ipt(conn, turn, ipt);

    connection_flush(conn); // Force flush
    LONGS_EQUAL(MSG_IPT_LEN, fakenetio_get_last_sent_size());

    fakeflukaio_get_last_sent_message(conn, &msg);
    MESSAGE_HEADER_EQUAL(&msg, N_IPT, sizeof(flukaio_ipt_data_t));
    LONGS_EQUAL(turn, msg.data.ipt.turn);
    LONGS_EQUAL(ipt, msg.data.ipt.ipt);

}

/* Send Particle */
TEST(FlukaIO, send_particle)
{
    n = flukaio_send_particle(conn, &expected);
    CHECK(n != -1);

    connection_flush(conn);// Force flush
    LONGS_EQUAL(MSG_PART_LEN, fakenetio_get_last_sent_size());

    fakeflukaio_get_last_sent_message(conn, &msg);
    MESSAGE_HEADER_EQUAL(&msg, N_PART, sizeof(expected));
    PARTICLES_EQUAL(&expected, &msg.data.particle);
}
/* Send Particle n times */
TEST(FlukaIO, send_particle_n_times)
{
    int times = 25;
    for(int i = 0; i < times; ++i) {
        expected.id = i;
        n = flukaio_send_particle(conn, &expected);
        CHECK(n != -1);
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


/* End Computation */
TEST(FlukaIO, send_end_comp)
{
    n = flukaio_send_eoc(conn);
    CHECK(n != -1);
    LONGS_EQUAL(MSG_EOC_LEN, fakenetio_get_last_sent_size());

    fakeflukaio_get_last_sent_message(conn, &msg);
    MESSAGE_HEADER_EQUAL(&msg, N_END, 0);
}

/* Send handshake */
TEST(FlukaIO, send_handshake)
{
    uint16_t major = 32;
    uint16_t minor = 55;
    uint32_t key = 45668;
    n = flukaio_send_hsk(conn, major, minor, key);
    LONGS_EQUAL(MSG_HSK_LEN, fakenetio_get_last_sent_size());

    fakeflukaio_get_last_sent_message(conn, &msg);
    MESSAGE_HEADER_EQUAL(&msg, N_HSK, sizeof(flukaio_hsk_data_t));
    LONGS_EQUAL(major, msg.data.hsk.major);
    LONGS_EQUAL(minor, msg.data.hsk.minor);
    LONGS_EQUAL(key,   msg.data.hsk.key);

}

/* Sync reads */
TEST(FlukaIO, wait_message_empty_buffer)
{
    n = flukaio_wait_message(conn, &msg);

    LONGS_EQUAL(-1, n);
}

TEST(FlukaIO, wait_with_empty_buffer_delay)
{

    // TODO: Add delay
    n = flukaio_wait_message(conn, &msg);

    LONGS_EQUAL(-1, n);
}

TEST(FlukaIO, wait_until_somethig_is_read)
{
    fakeflukaio_insert_particle_message(conn, &expected);

    n = flukaio_wait_message(conn, &msg);
    CHECK(n != -1);
}

TEST(FlukaIO, wait_message_particle)
{
    // wait particle message
    fakeflukaio_insert_particle_message(conn, &expected);

    n = flukaio_wait_message(conn, &msg);
    CHECK(n != -1);

    MESSAGE_HEADER_EQUAL(&msg, N_PART, sizeof(expected));
    PARTICLES_EQUAL(&expected, &msg.data.particle);
}


TEST(FlukaIO, wait_NULL_message_particle)
{
    // wait in NULL buffer skips next message
    fakeflukaio_insert_particle_message(conn, &expected);

    n = flukaio_wait_message(conn, NULL);
    CHECK(n != -1);
}

TEST(FlukaIO, wait_message_particle_twice)
{
    // Insert two messages and check they are read correctly
    fakeflukaio_insert_particle_message(conn, &expected);
    fakeflukaio_insert_particle_message(conn, &expected);

    n = flukaio_wait_message(conn, &msg);
    CHECK(n != -1);

    MESSAGE_HEADER_EQUAL(&msg, N_PART, sizeof(expected));
    PARTICLES_EQUAL(&expected, &msg.data.particle);

    memset(&msg, 0, sizeof(msg));
    n = flukaio_wait_message(conn, &msg);
    CHECK(n != -1);

    MESSAGE_HEADER_EQUAL(&msg, N_PART, sizeof(expected));
    PARTICLES_EQUAL(&expected, &msg.data.particle);
}

// Insert particle
TEST(FlukaIO, insert_particle)
{
    n = fakeflukaio_insert_particle_message(conn, &expected);

    LONGS_EQUAL(MSG_PART_LEN, n);
    LONGS_EQUAL(conn->in_buffer_len, n);
    fakeflukaio_clean_buffers(conn);
}

TEST(FlukaIO, insert_two_particles)
{
    n = fakeflukaio_insert_particle_message(conn, &expected);
    LONGS_EQUAL(MSG_PART_LEN, n);
    LONGS_EQUAL(conn->in_buffer_len, n);

    n = fakeflukaio_insert_particle_message(conn, &expected);
    LONGS_EQUAL(MSG_PART_LEN, n);
    LONGS_EQUAL(conn->in_buffer_len, n+n);

    fakeflukaio_clean_buffers(conn);
}
