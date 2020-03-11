import { app as _app, BrowserWindow as _BrowserWindow } from 'electron'
const app = _app
const BrowserWindow = _BrowserWindow

let mainWindow

function createWindow () {
  mainWindow = new BrowserWindow({width: 550, height: 650})
  mainWindow.loadURL('http://localhost:8000/index.html');
  mainWindow.on('closed', function () {
    mainWindow = null
  })
}

app.on('ready', createWindow)

app.on('window-all-closed', function () {
  app.quit()
});

app.on('activate', function () {
  if (mainWindow === null) {
    createWindow()
  }
})

function Supercharge() {
	var seq = document.getElementById("sequence-field").value;
	var site = document.getElementById("binding-field").value;
	var thres = document.getElementById("num-mu").value;

	console.log(typeof seq)
	
	eel.testing(seq,site,thres)(Set_result)
}

function Set_result(result) {
	document.getElementById("result-field").value = result;
}

var size = [window.width,window.height];

$(window).resize(function(){
    window.resizeTo(size[0],size[1]);
});