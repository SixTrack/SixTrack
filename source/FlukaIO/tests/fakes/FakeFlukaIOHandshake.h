#ifndef FLUKAIO_FAKEHANDSHAKE_H__
#define FLUKAIO_FAKEHANDSHAKE_H__

#include "FakeFlukaIO.h"

#ifdef __cplusplus
extern "C" {
#endif

extern int fakeflukaio_handshake_client_called;
extern int fakeflukaio_handshake_server_called;

int fakeflukaio_handshake_client(flukaio_connection_t * conn);
int fakeflukaio_handshake_server(flukaio_connection_t * conn);

#ifdef __cplusplus
}
#endif

#endif
