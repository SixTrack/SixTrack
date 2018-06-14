#ifndef FLUKAIO_FLUKAIO_H__
#define FLUKAIO_FLUKAIO_H__

#include <stdlib.h>
#include <stdint.h>

#include "Message.h"
#include "ParticleInfo.h"
#include "Connection.h"

#define FLUKAIO_MAJOR_VERSION 3
#define FLUKAIO_MINOR_VERSION 1
#define FLUKAIO_VERSION FLUKAIO_MAJOR_VERSION "." FLUKAIO_MINOR_VERSION

/** 
 * Maximum ratio of output buffer to be full before trying to send data
 * for example with 0.8 the output buffer will be sent when the buffer is full over the 
 * 80% of its capacity. 
 *
 * Set to 0 to disable caching
 */
#define FLUKAIO_OUT_CACHE 0.9
/** Maximum ratio of input buffer remaining before trying to read fron network.
 * For example with a value of 0.2 it will only try to read from network when the buffer is at
 * 20% or less of its capacity
 *
 * Set to 0 to disable caching
 */
#define FLUKAIO_IN_CACHE 0.1


#ifdef __cplusplus
extern "C" {
#endif

ssize_t flukaio_send_particle(flukaio_connection_t *conn, const particle_info_t *part);
ssize_t flukaio_send_ipt(flukaio_connection_t *conn, uint32_t turn, uint16_t ipt);
ssize_t flukaio_send_eob(flukaio_connection_t *conn);
ssize_t flukaio_send_eoc(flukaio_connection_t *conn);
ssize_t flukaio_send_int(flukaio_connection_t *conn, const int id, const int value);
ssize_t flukaio_send_double(flukaio_connection_t *conn, const int id, const double value);
ssize_t flukaio_send_npart(flukaio_connection_t *conn, const uint32_t npart);
ssize_t flukaio_send_hsk(flukaio_connection_t *conn, uint16_t major, uint16_t minor, uint32_t key);

ssize_t flukaio_receive_message(flukaio_connection_t *conn, flukaio_message_t *msg);
ssize_t flukaio_wait_message(flukaio_connection_t *conn, flukaio_message_t *msg);

flukaio_connection_t* flukaio_conn();
flukaio_connection_t* flukaio_connect(flukaio_connection_t *conn, const char *host, int port);
void flukaio_disconnect(flukaio_connection_t *conn);

extern int (*flukaio_handshake_client)(flukaio_connection_t * conn);

#ifdef __cplusplus
}
#endif

#endif
