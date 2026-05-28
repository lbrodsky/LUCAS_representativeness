#!/bin/bash -e

DB=lucas
DUMP_PATH=/data/db_st_lucas_dump.sql.7z

createdb -U postgres $DB
7z e -so $DUMP_PATH | sed '/^GRANT /d; /^ALTER DEFAULT PRIVILEGES /d' | psql -U postgres -d $DB

exit 0
