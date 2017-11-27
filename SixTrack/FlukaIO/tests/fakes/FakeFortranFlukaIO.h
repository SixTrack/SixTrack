#ifndef FLUKAIO_FAKEFORTRANFLUKAIO_H__
#define FLUKAIO_FAKEFORTRANFLUKAIO_H__

#include "Connection.h"
#include "FlukaIOServer.h"
#include "FortranFlukaIO.h"

#ifdef __cplusplus
extern "C" {
#endif

int store_connection(flukaio_connection_t *conn);
int store_server(flukaio_server_t *server);

#ifdef __cplusplus
}
#endif

#endif
