#ifndef FLUKAIO_FLUKAIO_PRIV_H__
#define FLUKAIO_FLUKAIO_PRIV_H__

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <errno.h>

#include "Connection.h"
#include "FlukaIO.h"

ssize_t flukaio_read_message_from(flukaio_message_t *msg, const void *buffer, size_t buffer_len);
ssize_t flukaio_write_pkt_bin(const char *data, const size_t datalen, char *buf, const size_t buflen);

#endif
