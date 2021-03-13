import cfdtools.api as api

def api_output_init():
    myapi = api.api_output()
    defmodes = myapi.get_modes()
    assert defmodes == myapi._default
    api1 = api.api_output('internal')
    api2 = api.api_output('internal')
    assert api1.get_modes()==api2.get_modes()
    ap1.set_default()
    assert defmodes == api1.get_modes()

