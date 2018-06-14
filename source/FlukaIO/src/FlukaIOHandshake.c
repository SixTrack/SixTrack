#include "FlukaIOHandshake.h"

int (*flukaio_handshake_client)(flukaio_connection_t * conn) = &impl_flukaio_handshake_client;
int (*flukaio_handshake_server)(flukaio_connection_t * conn) = &impl_flukaio_handshake_server;

/**
 * Handshake the server (used by the client)
 * checks if server implements a compatible version of the protocol
 * @param conn Connection object
 * @return true if versions are compatible, otherwise it disconnects
 */
int impl_flukaio_handshake_client(flukaio_connection_t * conn) {
	uint32_t key = 123141;
	ssize_t n = flukaio_send_hsk(conn, FLUKAIO_MAJOR_VERSION, FLUKAIO_MINOR_VERSION, key);
	if (n < 0) return -1;

	// Redundant check
	flukaio_message_t msg;
	n = flukaio_wait_message(conn, &msg);
	if (msg.type != N_HSK) return -1;
	if (msg.data.hsk.major != FLUKAIO_MAJOR_VERSION) return -1;
	if (msg.data.hsk.key   != key) return -1;
	return 0;
}

/**
 * Handshake the client (used by the server)
 * checks if client implements a compatible version of the protocol
 * @param conn Connection object
 * @return true if versions are compatible, otherwise it disconnects
 */
int impl_flukaio_handshake_server(flukaio_connection_t * conn) {
	flukaio_message_t msg;
	ssize_t n = flukaio_wait_message(conn, &msg);

	if (n < 0) return -1;
	if (msg.type != N_HSK) return -1;
	if (msg.data.hsk.major != FLUKAIO_MAJOR_VERSION) return -1;

	// Reply with our version plus the original key
	n = flukaio_send_hsk(conn, FLUKAIO_MAJOR_VERSION, FLUKAIO_MINOR_VERSION, msg.data.hsk.key);
	if (n < 0) return -1;
	return 0;
}

