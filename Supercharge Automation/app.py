import eel
eel.init('web')

@eel.expose
def my_python_function(param1, param2):
    print (param1 + param2)

eel.start('main.html', block=False)

eel.my_js_function('Hello ', 'World')

while True:
    eel.sleep(10)