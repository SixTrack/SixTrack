#ifndef FLUKAIO_CONNECTION_H__
#define FLUKAIO_CONNECTION_H__

#include <stdlib.h>
#include <sys/types.h>

#include "Message.h"

#ifdef __cplusplus
extern "C" {
#endif

#define IN_BUFFER_LEN 10240
#define OUT_BUFFER_LEN 10240

// These two values are specified in seconds
#define DEFAULT_READ_TIMEOUT 3000
#define DEFAULT_WRITE_TIMEOUT 3000

/** 
 * Connection information
 */
typedef struct _flukaio_connection_t {
	int fd;                             /**< Socket file descriptor */

	char   in_buffer[IN_BUFFER_LEN];    /**< Input data buffer */
	size_t in_buffer_size;              /**< Total size of the buffer (bytes) */
	size_t in_buffer_len;               /**< Length of the data in the buffer (bytes) */
	size_t in_buffer_start;             /**< Start position of the data in the buffer (bytes) */
	size_t in_buffer_end;               /**< End position of the data in the buffer (bytes) */

	char   out_buffer[OUT_BUFFER_LEN];  /**< Output data buffer */
	size_t out_buffer_size;             /**< Total size of the buffer (bytes) */
	size_t out_buffer_len;              /**< Length of the data in the buffer (bytes) */

	int read_timeout;                   /**< Timeout for reading operations (seconds) */
	int write_timeout;                  /**< Timeout for writing operations (seconds) */

	int     (*connect)(const char *host, int portnum);

	ssize_t (*read)(int fd, void *buf, size_t n);
	ssize_t (*write)(int fd, const void *buf, size_t len);

	int     (*can_write)(int fd, long timeout_sec, long timeout_usec);
	int     (*can_read)(int fd, long timeout_sec, long timeout_usec);

	int     (*set_nonblocking)(int fd);
	int     (*set_nodelay)(int fd);
} flukaio_connection_t;

flukaio_connection_t *connection_create(int fd);
int connection_connect(flukaio_connection_t *conn, const char *host, int port);

void    connection_destroy(flukaio_connection_t *conn);

ssize_t connection_push_message(flukaio_connection_t * conn,
		unsigned char type, const void *data, size_t datalen);
ssize_t connection_write(flukaio_connection_t *conn);
ssize_t connection_read(flukaio_connection_t *conn);
int     connection_can_read(flukaio_connection_t *conn);
int     connection_can_write(flukaio_connection_t *conn);

ssize_t connection_flush(flukaio_connection_t *conn);

int connection_set_read_timeout(flukaio_connection_t *conn, long timeout);
int connection_set_write_timeout(flukaio_connection_t *conn, long timeout);

ssize_t connection_write_message_bin(
		const unsigned char type,
		const void *data, const size_t datalen,
		void *buf, const size_t buflen);
ssize_t connection_receive_message(flukaio_connection_t *conn, flukaio_message_t *msg);
ssize_t connection_read_message_from(flukaio_message_t *msg, const void *buffer, size_t buffer_len);

#ifdef __cplusplus
}
#endif

#endif
