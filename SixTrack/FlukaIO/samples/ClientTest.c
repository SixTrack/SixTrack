#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <signal.h>
#include <string.h>
#include <assert.h>

#include <unistd.h>

#include "FlukaIO.h"

#define NPARTS 20000
#define NTURNS 150
#define PORT 14999

static flukaio_connection_t *conn = NULL;

void die(char *msg);
void sigint_handler(int sig_no);

int main(int argc, char *argv[]) {

	int port = PORT;
	long i, j;
	int n;
	flukaio_message_t msg;
	particle_info_t part;

	part.id =   0;
	part.gen =   1;
	part.weight = 1.0;
	part.x =   -6.319439221E-01;
	part.y =   1.147679085E+00;
	part.z =   1.0;
	part.tx =   -3.856529999E-05;
	part.ty =   -2.004572999E-06;
	part.tz =    5.004483739E-06;
	part.aa =   1;
	part.zz =   1;
	part.m =    450e0;
	part.pc =  4520e0;
	part.t =  5E-3;

	/* Install signal handler */
	struct sigaction action;
	memset(&action, 0, sizeof(action));
	action.sa_handler = &sigint_handler;
	action.sa_flags   = SA_RESTART;
	sigaction(SIGINT, &action, NULL);
	action.sa_handler = SIG_IGN;
	sigaction(SIGPIPE, &action, NULL);

	/* Main */
	conn = flukaio_connect(flukaio_conn(), "localhost", port);
	if (conn == NULL)
	{
		die("Error connecting to server, wrong protocol version?");
	}

	printf("Connected to server\n");
	for (j = 0; j < NTURNS; ++j)
	{
		printf("IPT %ld\n", j+1);
		n = flukaio_send_ipt(conn, j+1, 0);
		if (n < 0) {
			die("Error occurred sending IPT");
			break;
		}
		for (i = 0; i < NPARTS; ++i)
		{
			part.id = (unsigned int)i+1;
			n = flukaio_send_particle(conn, &part);
			if (n < 0) {
				die("Error occurred sending P");
				break;
			}
#ifdef VERBOSE
			printf(">");
#endif
			// Read incomming messages
			while (flukaio_receive_message(conn, &msg) > 0) {
#ifdef VERBOSE
				if (msg.type == N_PART)
					printf("<");
#endif
			}
		}
		n = flukaio_send_eob(conn);
		if (n < 0) {
			die("Error occurred sending EOB");
			break;
		}
#ifdef VERBOSE
		printf("EOB\n");
#endif
		// Read incomming messages until end of turn
		do {
			n = flukaio_wait_message(conn, &msg);
			if (n == -1) {
				die("Server timeout when waiting end of computation");
			}
#ifdef VERBOSE
			if (msg.type == N_PART)
				printf("<");
#endif
		} while (msg.type != N_EOB);
#ifdef VERBOSE
		printf("t\n");
#endif
		printf("Turn %ld finished\n", j+1);
	}

	printf("\nDone sending particles\n");
	printf("Closing connection...\n");

#ifdef VERBOSE
	printf("f\n");
#endif
	n = flukaio_send_eoc(conn);
	if (n < 0) {
		die("Error occurred sending");
	}

	n = flukaio_wait_message(conn, &msg);
	if (n == -1) {
		die("Server timeout when waiting end of computation");
	}
	if (msg.type != N_END) {
		die("Unexpected message received");
	}

	flukaio_disconnect(conn);
	conn = NULL;

	return 0;
}

/* Close */
void die(char *msg) {
	perror(msg);

	if (conn)
	{
		flukaio_disconnect(conn);
		conn = NULL;
	}

	exit(1);
}

/* Signal handlers */
void sigint_handler(int sig_no)
{
	die("Interrupted: Ending connection");
}

