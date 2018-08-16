#ifndef FLUKAIO_FAKEFLUKAIOSERVER_H__
#define FLUKAIO_FAKEFLUKAIOSERVER_H__

#include "FlukaIOServer.h"
#include "FakeNetIO.h"
#include "FakeConnection.h"

#ifdef __cplusplus
extern "C" {
#endif

flukaio_server_t* fakeflukaio_server_create();

#ifdef __cplusplus
}
#endif

#endif
