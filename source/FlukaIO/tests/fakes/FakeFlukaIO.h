#ifndef FLUKAIO_FAKEFLUKAIO_H__
#define FLUKAIO_FAKEFLUKAIO_H__

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stddef.h>

#include "FlukaIO.h"
#include "Connection.h"

#ifdef __cplusplus
extern "C" {
#endif

int flukaio_send_message(flukaio_connection_t *conn, const flukaio_message_t *msg);

ssize_t fakeflukaio_insert_message(
        flukaio_connection_t *conn,
        const unsigned char type,
        const void *data, const size_t datalen);
ssize_t fakeflukaio_insert_eob_message(flukaio_connection_t *conn);
ssize_t fakeflukaio_insert_ipt_message(flukaio_connection_t *conn,
		const uint32_t turn, const uint16_t ipt);
ssize_t fakeflukaio_insert_particle_message(flukaio_connection_t *conn,
		const particle_info_t *part);
ssize_t fakeflukaio_insert_eoc_message(flukaio_connection_t *conn);
ssize_t fakeflukaio_insert_config_message(flukaio_connection_t *conn);
ssize_t fakeflukaio_insert_unknown_message(flukaio_connection_t *conn);
ssize_t fakeflukaio_insert_hsk_message(flukaio_connection_t *conn,
		const uint16_t major, const uint16_t minor, const uint32_t key);

ssize_t fakeflukaio_get_last_sent_message(flukaio_connection_t *conn, flukaio_message_t *msg);

void fakeflukaio_clean_buffers(flukaio_connection_t *conn);

// Defined as private in FlukaIO_private.h
ssize_t connection_read_message_from(flukaio_message_t *msg, const void *buffer, size_t buffer_len);
ssize_t connection_write_message_bin(
        const flukaio_message_type_t type,
        const void *data, const size_t datalen,
        void *buf, const size_t buflen);

ssize_t flukaio_flush(flukaio_connection_t *conn);

#ifdef __cplusplus
}
#endif

#endif
