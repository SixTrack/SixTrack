#include <stdlib.h>
#include <stdio.h>

#include <unistd.h>
#include <string.h>
#include <signal.h>

#include "FlukaIO.h"
#include "FlukaIOServer.h"

#define PORT 14999

static flukaio_connection_t *conn = NULL;
static flukaio_server_t *server = NULL;
void die(char *msg);
void sigint_handler(int sig_no);

/* */
int main(int argc, char *argv[]) {

	int n;
	flukaio_message_t msg;
	uint32_t turn  = 0;
	uint32_t count = 0;
	uint32_t total = 0;

	/* Install signal handlers */
	struct sigaction action;
	memset(&action, 0, sizeof(action));
	action.sa_handler = &sigint_handler;
	action.sa_flags   = SA_RESTART;
	sigaction(SIGINT, &action, NULL);
	action.sa_handler = SIG_IGN;
	sigaction(SIGPIPE, &action, NULL);

	server = flukaio_server_create();
	if (!server)
	{
		die("** Error, create server object\n");
	}

	n = flukaio_server_start(server, PORT);
	if (n == -1)
	{
		die("** Error, could not start server\n");
	}

	printf("Listening on port %d\n", n);

	while(1) {

		printf("Waiting new connection...\n");
		conn = flukaio_server_accept(server);
		if (!conn) {
			printf("Failed connection attempt, wrong protocol version?\n");
		} else {
			printf("New connection accepted\n");
			count = 0;
			//connection_set_read_timeout(conn, 20);

			while (1) {
				n = flukaio_wait_message(conn, &msg);

				if (n < 0) {
					printf("** Client timeout\n");
					break;
				}

				if (msg.type == N_IPT) {
					turn = msg.data.ipt.turn;
					printf("Starting turn %d\n", turn);
				}
				else if (msg.type == N_PART) {
					count++;
					n = flukaio_send_particle(conn, &msg.data.particle);
					if (n < 0) {
						printf("** Error occurred sending\n");
						break;
					}
				}
				else if (msg.type == N_EOB) {
					printf("End of turn %d, received %d particles\n", turn, count);
					flukaio_send_eob(conn);
					total += count;
					count = 0;
				}
				else if (msg.type == N_END) {
					flukaio_send_eoc(conn);
					printf("Connection finished\n");
					break;
				}
			}

			flukaio_disconnect(conn);
			conn = NULL;
		}
	}

	flukaio_server_shutdown(server);

	return 0;
}

/* Close */
void die(char *msg) {
	printf(msg);

	if(conn)
	{
		flukaio_disconnect(conn);
		conn = NULL;
	}
	flukaio_server_shutdown(server);

	exit(1);
}

/* Signal handlers */
void sigint_handler(int sig_no)
{
	die("\n** Shuting down server...");
}

