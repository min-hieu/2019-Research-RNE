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