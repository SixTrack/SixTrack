#include "Connection.h"
#include "fakes/FakeNetIO.h"
#include "fakes/FakeConnection.h"

#include "CommonTest.h"

TEST_GROUP(Connection)
{
    static const int fd = -1;
    flukaio_connection_t *conn;
    ssize_t oldretval;

    void setup()
    {
        oldretval = connection_read(NULL);

        conn = fakeconnection_create(fd);

        if (conn == NULL)
            FAIL("Connection should not be NULL now");
    }

    void teardown()
    {
        connection_destroy(conn);
        conn = NULL;
        fakenetio_set_read_retval(oldretval);
    }

};

TEST(Connection, ConnectionCreate)
{
    CHECK(conn != NULL);
    LONGS_EQUAL(fd, conn->fd);
    CHECK(conn->in_buffer != NULL);
    CHECK(conn->out_buffer != NULL);
    LONGS_EQUAL(0, conn->in_buffer_len);
    LONGS_EQUAL(0, conn->out_buffer_len);
    CHECK(conn->in_buffer_size != 0);
    CHECK(conn->out_buffer_size != 0);
    LONGS_EQUAL(DEFAULT_READ_TIMEOUT, conn->read_timeout);
    LONGS_EQUAL(DEFAULT_WRITE_TIMEOUT, conn->write_timeout);
}

TEST(Connection, connection_read_update_buffer_len)
{
    size_t prevlen = conn->in_buffer_len;
    fakenetio_set_read_retval(1);
    connection_read(conn);
    LONGS_EQUAL(prevlen + 1, conn->in_buffer_len);

}

TEST(Connection, connection_read_doesnt_update_buffer_len_if_error)
{
    size_t prevlen = conn->in_buffer_len;
    fakenetio_set_read_retval(-1);
    connection_read(conn);
    LONGS_EQUAL(prevlen, conn->in_buffer_len);

}

TEST(Connection, connection_set_read_timeout)
{
    connection_set_read_timeout(conn, 137);
    LONGS_EQUAL(conn->read_timeout, 137);
}

// Write message bin
//
TEST(Connection, write_message_bin_just_buffer)
{
	uint8_t data[8] = { 230, 123, 234, 54, 1, 2, 3, 4 };
	uint8_t buffer[MSG_HEADER_LEN+sizeof(data)];
	uint8_t type = 2;
    int n = connection_write_message_bin(type, data, sizeof(data), buffer, sizeof(buffer));

    LONGS_EQUAL(MSG_HEADER_LEN+8, n);
}
TEST(Connection, write_message_bin_big_buffer)
{
	uint8_t data[8] = { 230, 123, 234, 54, 1, 2, 3, 4 };
	uint8_t buffer[1024];
	uint8_t type = 2;
    int n = connection_write_message_bin(type, data, sizeof(data), buffer, sizeof(buffer));

    LONGS_EQUAL(MSG_HEADER_LEN+8, n);
}
TEST(Connection, write_message_bin_small_buffer)
{
	uint8_t data[8] = { 230, 123, 234, 54, 1, 2, 3, 4 };
	uint8_t buffer[MSG_HEADER_LEN];
	uint8_t type = 2;
    int n = connection_write_message_bin(type, data, sizeof(data), buffer, sizeof(buffer));

    LONGS_EQUAL(-1, n);
}

