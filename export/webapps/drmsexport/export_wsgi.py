# can test out this wsgi entry point by providing arguments to uwsgi:
#   uwsgi --socket 171.64.103.154:50000 --protocol=http -w export_wsgi:export
# ultimately, we will create and provide to uwsgi an .ini file that contains the desired parameters

from drms_parameters import DRMSParams

# import the flask app from export.py
from export import export, APP_DEBUG, APP_LOG

if __name__ == "__main__":
    drms_params = DRMSParams()
    export.run(host=drms_params.EXPORT_WEB_SERVER, port=drms_parms.EXPORT_WEB_SERVER_PORT)
else:
    pass
