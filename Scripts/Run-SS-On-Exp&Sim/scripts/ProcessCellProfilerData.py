class Dict2Class(object):
    """
        goal:  Change Dictionary Into Class
               this function is useful to read experimental output pickle files.

    """

    def __init__(self, my_dict):
        for key in my_dict:
            setattr(self, key, my_dict[key])