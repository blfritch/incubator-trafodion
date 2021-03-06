//------------------------------------------------------------------
//
// @@@ START COPYRIGHT @@@
//
// (C) Copyright 2006-2014 Hewlett-Packard Development Company, L.P.
//
//  Licensed under the Apache License, Version 2.0 (the "License");
//  you may not use this file except in compliance with the License.
//  You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
//
// @@@ END COPYRIGHT @@@

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "seabed/fserr.h"
#include "seabed/ms.h"
#include "seabed/pctl.h"
#include "seabed/pevents.h"

#include "tchkfe.h"
#include "tmsfsutil.h"
#include "tutil.h"


int main(int argc, char *argv[]) {
    int     exit;
    int     ferr;
    int     lerr;
    char    name[BUFSIZ];
    char   *p;
    char    recv_buffer[BUFSIZ];
    MS_SRE  sre;

    msfs_util_init_role(false, &argc, &argv, msg_debug_hook);
    ferr = msg_mon_process_startup(true); // system messages
    TEST_CHK_FEOK(ferr);
    msg_mon_enable_mon_messages(true);
    ferr = msg_mon_get_my_process_name(name, sizeof(name));
    TEST_CHK_FEOK(ferr);
    exit = 0;
    p = fgets(recv_buffer, BUFSIZ, stdin);
    assert(p != NULL);
    printf(recv_buffer);
    while (!exit) {
        do {
            lerr = XWAIT(LREQ, -1);
            TEST_CHK_WAITIGNORE(lerr);
            lerr = XMSG_LISTEN_((short *) &sre, // sre
                                0,              // listenopts
                                0);             // listenertag
        } while (lerr == XSRETYPE_NOWORK);
        ferr = XMSG_READDATA_(sre.sre_msgId,  // msgid
                              recv_buffer,    // reqdata
                              BUFSIZ);        // bytecount
        util_check("XMSG_READDATA_", ferr);
        if (sre.sre_flags & XSRE_MON) {
            MS_Mon_Msg *msg = (MS_Mon_Msg *) recv_buffer;
            if (msg->type == MS_MsgType_Close)
                exit = 1;
        } else {
            printf("srv %s: NON-mon message\n", name);
        }
        XMSG_REPLY_(sre.sre_msgId,       // msgid
                    NULL,                // replyctrl
                    0,                   // replyctrlsize
                    NULL,                // replydata
                    0,                   // replydatasize
                    0,                   // errorclass
                    NULL);               // newphandle
    }
    printf("server %s shutting down\n", name);
    ferr = msg_mon_process_shutdown();
    TEST_CHK_FEOK(ferr);
    return 0;
}
