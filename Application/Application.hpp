
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <cstdint>
#include <random>
#include <unistd.h>

#include "../CONFIG.h"
#include "Packet.hpp"
#include "../Untrusted/OLib.hpp"
#include "../Untrusted/PathOHeap.hpp"
#include "../Untrusted/SN.hpp"
#include "../Untrusted/RS.hpp"
#include "../Untrusted/BOS.hpp"
#include "../Untrusted/TC.hpp"
#include "../Untrusted/OS_RSandNS.hpp"
#include "../Untrusted/WN.hpp"
#include "../Untrusted/BN.hpp"

#include <openssl/ec.h>
#include <openssl/ecdh.h>
#include <openssl/ecdsa.h>
#include <openssl/conf.h>
#include <openssl/evp.h>
#include <openssl/err.h>
#include <openssl/obj_mac.h>

#define CLOCKS_PER_MS (CLOCKS_PER_SEC/1000)
