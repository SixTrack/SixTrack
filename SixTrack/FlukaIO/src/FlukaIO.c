#include "FlukaIO_private.h"

flukaio_connection_t* flukaio_conn() {
	return connection_create(-1);
} // flukaio_conn

/**
 * Blocks until a message is read.
 * Stores read message in msg
 * @return number bytes read, -1 if error (timeout, disconnect, etc.)
 */
ssize_t flukaio_wait_message(flukaio_connection_t *conn, flukaio_message_t *msg) {
	ssize_t n;
	int retries = 0;

	if (conn == NULL)
		return -1;

	while (retries < 2 && (n = flukaio_receive_message(conn, msg)) == -1) {
		// blocks waiting for data
		int retval = connection_can_read(conn);

		if (retval <= 0)
			return -1;

		retries++;
	}
	return n;
} // flukaio_wait_message

/**
 * Returns a message if inmediately available
 * Stores read message in msg
 * @return bytes read, -1 if error (disconnect, etc.)
 */
ssize_t flukaio_receive_message(flukaio_connection_t *conn, flukaio_message_t *msg)
{
	if (conn == NULL)
		return -1;

	if (conn->in_buffer_len <= conn->in_buffer_size * FLUKAIO_IN_CACHE)
		connection_read(conn);

	int n = connection_receive_message(conn, msg);
	if (n > 0) {
		if (msg) {
			switch(msg->type) {
				case N_PART:
				case N_EOB:
				case N_IPT:
				case N_END:
				case N_CONF:
				case N_HSK:
				case N_INT:
			        case N_DBLE:
					break;
				default:
					msg->type = N_ERR;
					return -1;
					break;
			}
		}
		return n;
	}

	if (msg) msg->type = N_ERR;
	return -1;
} // flukaio_receive_message

/**
 * Send particle
 * @param part particle to send
 * @return bytes sent, -1 if error
 */
ssize_t flukaio_send_particle(flukaio_connection_t *conn, const particle_info_t *part)
{
	ssize_t n = connection_push_message(conn, N_PART, part, sizeof(particle_info_t));

	if (n < 0) return n;

	if(conn->out_buffer_len >= conn->out_buffer_size * FLUKAIO_OUT_CACHE)
		connection_flush(conn);

	return n;

} // flukaio_send_particle

/**
 * Send End Of Batch
 * @return bytes sent, -1 if error
 */
ssize_t flukaio_send_eob(flukaio_connection_t *conn)
{
	ssize_t n = connection_push_message(conn, N_EOB, NULL, 0);
	connection_flush(conn);
	return n;
} // flukaio_send_eob

/**
 * Send Insertion Point
 * @param turn Current turn number
 * @param ipt Insertion point number
 * @return bytes sent, -1 if error
 */
ssize_t flukaio_send_ipt(flukaio_connection_t *conn, uint32_t turn, uint16_t ipt)
{
	flukaio_ipt_data_t data;
	data.turn = turn;
	data.ipt = ipt;

	ssize_t n = connection_push_message(conn, N_IPT, &data, sizeof(data));
	return n;
} // flukaio_send_ipt

/**
 * Send Handshake
 */
ssize_t flukaio_send_hsk(flukaio_connection_t *conn, uint16_t major, uint16_t minor, uint32_t key)
{
	flukaio_hsk_data_t data;
	data.minor = minor;
	data.major = major;
	data.key = key;

	ssize_t n = connection_push_message(conn, N_HSK, &data, sizeof(data));
	connection_flush(conn);
	return n;
} // flukaio_send_hsk

/**
 * Send End Of Computation
 * @return bytes sent, -1 if error
 */
ssize_t flukaio_send_eoc(flukaio_connection_t *conn)
{
	ssize_t n = connection_push_message(conn, N_END, NULL, 0);
	connection_flush(conn);
	return n;
} // flukaio_send_eoc

/**
 * Send an int message
 * @param conn   connection structrure
 * @param id     message id
 * @param value  integer value to send
 * @return bytes sent, -1 if error
 */
ssize_t flukaio_send_int(flukaio_connection_t *conn, const int id, const int value)
{
	flukaio_int_data_t data;

	data.id    = id;
	data.value = value;

	ssize_t n = connection_push_message(conn, N_INT, &data, sizeof(data));
	connection_flush(conn);
	return n;
} // flukaio_send_int

/**
 * Send a double message
 * @param conn   connection structrure
 * @param id     message id
 * @param value  double value to send
 * @return bytes sent, -1 if error
 */
ssize_t flukaio_send_double(flukaio_connection_t *conn, const int id, const double value)
{
	flukaio_dble_data_t data;

	data.id    = id;
	data.value = value;

	ssize_t n = connection_push_message(conn, N_DBLE, &data, sizeof(data));
	connection_flush(conn);
	return n;
} // flukaio_send_double

/**
 * Send maximum Number of particles
 * @return bytes sent, -1 if error
 */
ssize_t flukaio_send_npart(flukaio_connection_t *conn, const uint32_t npart)
{
	return flukaio_send_int(conn, N_NPART, npart);
} // flukaio_send_npart

/**
 * Helper routine, connect to server at host and port
 * @return 0 if success, -1 if error
 */
flukaio_connection_t *flukaio_connect(flukaio_connection_t *conn, const char *host, int port) {
	if (connection_connect(conn, host, port) < 0) return conn;

	if (flukaio_handshake_client(conn) < 0) {
		// Handshake failed
		errno = EPROTO;
		flukaio_disconnect(conn);
		return NULL;
	}

	return conn;
} // flukaio_connect

/**
 * Disconnect from server
 * @param conn Connection object
 */
void flukaio_disconnect(flukaio_connection_t *conn) {
	connection_destroy(conn);
} // flukaio_disconnect
