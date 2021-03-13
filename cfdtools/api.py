

class api_output():
    """class to handle library outputs
    """
    _available = ['internal', 'error', 'warning', 'standart', 'debug']
    _default = ['internal', 'error', 'warning', 'standart' ]

    def __init__(self, list=_default):
        self._api_output = []
        self.set_modes(list)

    def set_modes(self, list):
        list = [].append(list) # ensure list if arg is only an item
        check = all(mode in self._available for mode in list)
        if check:
            self._api_output = list
        else:
            self.print('internal','some output modes are unknown')
        return check

    def get_modes(self):
        return self._api_output

    def set_default(self):
        self.set_modes(self.default)

    def print(self, mode, *args):
        if mode in self._api_output:
            print(mode+': ',*args)