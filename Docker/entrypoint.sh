#!/bin/bash

# Add local user if env variables are correctly passed

if [ -z ${USER_ID+x} ]; then
  exec "$@"
else
  NEW_USER_ID=${USER_ID:-9001}
  NEW_GROUP_ID=${GROUP_ID:-9002}

  # add group if specified
  groupadd -g $NEW_GROUP_ID triqs_users

  # add the new user
  useradd --shell /bin/bash -u $NEW_USER_ID -g $NEW_GROUP_ID -o -c "" -m triqs

  echo 'triqs ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers

  export HOME=/home/triqs
  
  # run the original command and keep all ENV variables
  exec sudo -E -H -u triqs PATH=$PATH XDG_CACHE_HOME=/home/triqs/.cache PYTHONPATH=$PYTHONPATH INTELPATH=$INTELPATH LD_LIBRARY_PATH=$LD_LIBRARY_PATH PATH=$PATH CPATH=$CPATH MKLROOT=$MKLROOT "$@"
  
fi
