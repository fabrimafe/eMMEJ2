"""
This module is a helper module that contains functions that I 
use to while developing other modules (such as type chacking)

"""

def type_checker(arg_dict: dict, input_arg: list, mode: str):
    """
    Checks that the arguments are of the right type, if not,
    raise a TypeError

    call it as follow:
    if inside a class:
        type_chacker(self.__dict__, self.__init__.__annotations__.values())
    if calling on a function:
        type_chacker(some_func.__annotations__, locals().values())
    Args:
        arg_dict (dict): a dictionary containing all the inputs and their keys
        input_arg (list): only the values (in this case, the typing annotations)
            of the  __annotations__ (if in function: function.__annotations__.values()
                                    if in class: self.function.__annotations__.values())

    """
    # print('type_chacker activated')
    if (mode=='class'):
        for arg_key in arg_dict.keys():
            # print(f'arg_key: {arg_key}, arg_val: {arg_val}, inp_key: {inp_key}, inp_val: {inp_val}')
            # print(arg_dict[inp_key])
            # print((type(arg_dict[inp])))
            print(arg_key)
            print(input_arg[arg_key])
            print(type(arg_dict[arg_key]))
            
            if (type(arg_dict[arg_key]) != input_arg[arg_key]):
                raise TypeError(f"""{arg_key} if of wrong type
                expected: {input_arg[arg_key]}, got: {type(arg_dict[arg_key])}, \n   INPUT: {arg_dict[arg_key]}""")

    if (mode=='func'):
        for key,val in zip(arg_dict, input_arg):
            if (val != type(input_arg[key])):
                raise TypeError(f"""{key} if of wrong type
                expected: {val}, got: {type(input_arg[key])}, \n    INPUT: {val}, {type(val)}""")
        