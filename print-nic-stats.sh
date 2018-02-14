#! /bin/bash

ethtool -S eth2 | grep rx_dropped
ethtool -S eth2 | grep vport_dropped
ethtool -S eth2 | grep rx_fifo_errors
ethtool -S eth2 | grep rx_missed_errors
