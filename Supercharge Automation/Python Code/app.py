import eel
eel.init('web')

eel.start('index.html', size=(800,200), block=False)

eel.my_js_function('Hello ', 'World')

while True:
    eel.sleep(10)