#ifndef FLUKAIO_FAKECONNECTION_H__
#define FLUKAIO_FAKECONNECTION_H__

#include "Connection.h"

#ifdef __cplusplus
extern "C" {
#endif

flukaio_connection_t* fakeconnection_create(int fd);

#ifdef __cplusplus
}
#endif

#endif
