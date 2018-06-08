#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <assert.h>

#include "Connection.h"
#include "NetIO.h"


/**
 * Connect to a host and port
 * @return -1 if error
 */
int connection_connect(flukaio_connection_t *conn, const char *host, int port) {
	conn->fd = conn->connect(host, port);
	return conn->fd;
}

/**
 * Creates a connection object given a socket (file descriptor)
 * @return returns a connection, NULL if memory error
 */
flukaio_connection_t *connection_create(int fd)
{
	flukaio_connection_t *conn    = (flukaio_connection_t *)malloc(sizeof(flukaio_connection_t));
	if (conn)
	{
		conn->fd                    = fd;

		conn->in_buffer_size        = sizeof(conn->in_buffer);
		conn->out_buffer_size       = sizeof(conn->out_buffer);
		conn->in_buffer_len         = 0;
		conn->in_buffer_start       = 0;
		conn->in_buffer_end         = 0;
		conn->out_buffer_len        = 0;
		conn->read_timeout          = DEFAULT_READ_TIMEOUT;
		conn->write_timeout         = DEFAULT_WRITE_TIMEOUT;

		conn->connect = &netio_connect;
		conn->read = &netio_read;
		conn->write = &netio_write;
		conn->can_read = &netio_can_read;
		conn->can_write = &netio_can_write;
		conn->set_nonblocking = &netio_set_nonblocking;
		conn->set_nodelay = &netio_set_nodelay;
	}

	return conn;
}

/**
 * Closes connection and frees structures
 */
void connection_destroy(flukaio_connection_t *conn)
{
	if (conn == NULL)
		return;

	if(conn->fd != -1)
	{
		shutdown(conn->fd, SHUT_RDWR);
		close(conn->fd);
	}
	free(conn);
}

/**
 * Send message through connection
 * @return bytes sent, -1 if error
 */
ssize_t connection_push_message(flukaio_connection_t * conn,
		unsigned char type, const void *data, size_t datalen)
{
	if (conn == NULL)
		return -1;

	size_t out_buffer_space = conn->out_buffer_size - conn->out_buffer_len;
	ssize_t n = connection_write_message_bin(
			type, data, datalen,
			conn->out_buffer + conn->out_buffer_len,
			out_buffer_space);

	if (n <= 0) return n;
	conn->out_buffer_len += n;

	return n;
}

/**
 * Write message in buffer
 * @return number of bytes written, -1 if error
 */
ssize_t connection_write_message_bin(
		const unsigned char type,
		const void *data, const size_t datalen,
		void *buf, const size_t buflen)
{
	if (buf == NULL) return -2;
	if (data == NULL && datalen > 0) return -3;
	if (buflen < (datalen + MSG_HEADER_LEN)) return -4;

	void *bufp = buf;
	flukaio_message_t msg;

	msg.size = MSG_HEADER_LEN + datalen;
	msg.type = type;

	memcpy(bufp, &msg, MSG_HEADER_LEN);
	memcpy(bufp+MSG_HEADER_LEN, data, datalen);

	return msg.size;
}

/**
 * Returns a message if inmediately available
 * Stores read message in msg
 * @return bytes read, -1 if error (disconnect, etc.)
 */
ssize_t connection_receive_message(flukaio_connection_t *conn, flukaio_message_t *msg)
{
	if (conn == NULL)
		return -1;

	ssize_t len = connection_read_message_from(msg, conn->in_buffer + conn->in_buffer_start, conn->in_buffer_len);

	if (len > 0)
	{
		conn->in_buffer_len -= len;
		conn->in_buffer_start += len;
	}

	return len;
}


/**
 * Write/transmit output buffer, blocks until done
 * @return number of bytes written, -1 if connection error (disconnect, etc.)
 */
ssize_t connection_write(flukaio_connection_t *conn)
{
	if (conn == NULL)
		return -1;

	ssize_t n;
	n = conn->write(conn->fd, conn->out_buffer, conn->out_buffer_len);
	// write always sends the full buffer or reports an error
	if (n > 0)
		conn->out_buffer_len -= n;
	return n;
}

/**
 * Write all contents of connection
 * exactly the same as connection_write
 */
ssize_t connection_flush(flukaio_connection_t *conn)
{
	return connection_write(conn);
}

/**
 * Read from network and store in input buffer
 * @return number of read bytes, -1 if connection error (disconnect, etc.)
 */
ssize_t connection_read(flukaio_connection_t *conn)
{
	if (conn == NULL)
		return -1;

	if (conn->in_buffer_start > (conn->in_buffer_size>>1)) {
		memmove(conn->in_buffer, conn->in_buffer + conn->in_buffer_start, conn->in_buffer_len);
		conn->in_buffer_start = 0;
		conn->in_buffer_end = conn->in_buffer_len;
	}
	ssize_t n = conn->read(conn->fd,
		conn->in_buffer + conn->in_buffer_end,
		conn->in_buffer_size-conn->in_buffer_end);
	if (n > 0) {
		conn->in_buffer_end += n;
		conn->in_buffer_len += n;
	}
	return n;
}

/**
 * @return 1 if can read, -1 on error
 */
int connection_can_read(flukaio_connection_t *conn)
{
	if (conn == NULL)
		return -1;

	return conn->can_read(conn->fd, conn->read_timeout, 0);
}

/**
 * @return 1 if can write, -1 on error
 */
int connection_can_write(flukaio_connection_t *conn)
{
	if (conn == NULL)
		return -1;

	return conn->can_write(conn->fd, conn->write_timeout, 0);
}

/**
 * Sets the read timeout for this specific connection
 * @param time to wait in seconds
 * @return set timeout, -1 if error
 */
int connection_set_read_timeout(flukaio_connection_t *conn, long timeout)
{
	if (conn == NULL)
		return -1;

	return conn->read_timeout = timeout;
}

int connection_set_write_timeout(flukaio_connection_t *conn, long timeout)
{
	if (conn == NULL)
		return -1;

	return conn->write_timeout = timeout;
}

/**
 * read message from a given buffer
 * Stores read message in msg
 * @return bytes read, -1 if error (disconnect, etc.)
 */
inline ssize_t connection_read_message_from(flukaio_message_t *msg, const void *buffer, size_t buffer_len)
{
	flukaio_message_t *mymsg = msg;
	flukaio_message_t m;

	if (mymsg == NULL)
		mymsg = &m;

	if (buffer == NULL)
		return -1;

	if (buffer_len < MSG_HEADER_LEN)
		return -1;

	memcpy(mymsg, buffer, MSG_HEADER_LEN);

	// Full data not ready
	if (buffer_len < mymsg->size)
		return -1;

	if (msg != NULL)
	{
		memcpy(msg, buffer, mymsg->size);
	}

	return mymsg->size;

}

