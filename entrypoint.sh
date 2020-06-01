#!/usr/bin/env bash

set -e

# compile the frontend app
npm run-script build --prefix /code/src/genui_gui/

# migrate the database and set everything up
python manage.py migrate --noinput
python manage.py genuisetup
python manage.py collectstatic --no-input

# ensure that only the genui user group has access to media files
chgrp -R ${GENUI_USER_GROUP} /code/src/genui/media
chmod g+rwxs,a+xs /code/src/genui/media # FIXME: here we allow all users to list directories and read the media files -> less permission headaches setting up the web server, but it is a security concern in some deployment settings
umask 002

# allow the genui user group to bind the 443 port
touch /etc/authbind/byport/443
chgrp ${GENUI_USER_GROUP} /etc/authbind/byport/443
chmod g+rwx /etc/authbind/byport/443

# make sure the user has a home directory in the container and has the required permissions
mkdir -p /home/${GENUI_USER}
chgrp -R ${GENUI_USER_GROUP} /home/${GENUI_USER}
chown -R ${GENUI_USER} /home/${GENUI_USER}
chmod 770 /home/${GENUI_USER}

# execute the command
exec runuser -u  ${GENUI_USER} -- "$@"