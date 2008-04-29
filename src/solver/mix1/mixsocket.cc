/*****************************************************************************/
/*              8888888         88888888         88888888                    */
/*            8                8                8                            */
/*           8                 8                8                            */
/*           8                  88888888         88888888                    */
/*           8      8888                8                8                   */
/*            8       8                 8                8                   */
/*              888888         888888888        888888888                    */
/*                                                                           */
/*       A Two-Dimensional General Purpose Semiconductor Simulator.          */
/*                                                                           */
/*  GSS 0.4x                                                                 */
/*  Last update: March 29, 2007                                              */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/
#include "mix1.h"
#include "log.h"
#include "string.h"


sigjmp_buf net_buf;

/* This function gets called whenever the pipe broken.
See the signal(2) manpage for more information. */
inline void signal_handler(int signum)
{
  switch (signum)
  {
  case SIGPIPE:
    printf("Received signal SIGPIPE. It seems ngspice terminated.\n");
    siglongjmp(net_buf,1);
  default:
    break;
  }
}

/* ----------------------------------------------------------------------------
 * NDEVServer:  This function set TCP/IP socket and connect to NGSPICE. 
 */
int NDEVServer(int port,int &listener, int &client)
{
  char dotted_ip[15]; /* Buffer for converting the resolved address to a readable format. */
  struct sockaddr_in sa; /* Connection address. */
  socklen_t sa_len;
  char sendbuf[128];
  signal(SIGPIPE, signal_handler);
  listener = socket(PF_INET, SOCK_STREAM, IPPROTO_IP);
  if (listener < 0)
  {
    gss_log.string_buf()<<"Unable to create a listener socket: "<<strerror(errno)<<"\n";
    gss_log.record();
    return 1;
  }
  /* Now bind the listener to a local address. This uses
   the same sockaddr_in structure as connect. */
  sa_len = sizeof(sa);
  //bzero(&sa, sa_len);
  memset(&sa,0,sa_len);
  sa.sin_family = AF_INET;
  sa.sin_port = htons(port);
  sa.sin_addr.s_addr = htonl(INADDR_ANY); /* Listen on all interfaces. */
  if (bind(listener, (sockaddr *)&sa, sa_len) < 0)
  {
    gss_log.string_buf()<<"Unable to bind to port "<<port<<": "<<strerror(errno)<<"\n";
    gss_log.record();
    return 1;
  }
  /* Let the networking system know we're accepting
  connections on this socket. Ask for a connection
  queue of five clients. (If more than five clients
  try to connect before we call accept, some will
  be denied.) */
  if (listen(listener, 1) < 0)
  {
    gss_log.string_buf()<<"Unable to listen: "<<strerror(errno)<<"\n";
    gss_log.record();
    return 1;
  }
  gss_log.string_buf()<<"Waiting for ngspice...";
  gss_log.record();
  client = accept(listener, (sockaddr *)&sa, &sa_len);
  /* We now have a live client. Print information
    about it and then send something over the wire. */
  inet_ntop(AF_INET, &sa.sin_addr, dotted_ip, 15);

  /* Waiting for ngspice query data */
  recv(client,sendbuf,128,MSG_FLAG);
  if(strcmp(sendbuf,NG_QUERY)) /* Query error */
  {
    gss_log.string_buf()<<"Query information error from ngspice.\n";
    gss_log.record();
    return 1;
  }
  gss_log.string_buf()<<"Received ngspice connection from "<<dotted_ip<<"\n";
  gss_log.record();
  /* ACK to ngspice */
  sprintf(sendbuf,NDEV_REPLY);
  send(client,sendbuf,128,0);

  return 0;
}
