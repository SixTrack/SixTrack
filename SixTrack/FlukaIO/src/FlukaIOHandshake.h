#ifndef FLUKAIO_HANDSHAKE_H__
#define FLUKAIO_HANDSHAKE_H__

#include "FlukaIO.h"

#ifdef __cplusplus
extern "C" {
#endif

int impl_flukaio_handshake_client(flukaio_connection_t * conn);
int impl_flukaio_handshake_server(flukaio_connection_t * conn);

#ifdef __cplusplus
}
#endif

#endif
