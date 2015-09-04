var lastHash = window.location.hash;

function showAdvanced(element) {
	if ($(".advancedLeft").is(":visible")) {
		$(".advancedLeft, .advancedRight").slideUp(400);	
		$(element).html("Show advanced settings");
	} else {
		$(".advancedLeft, .advancedRight").slideDown(400);	
		$(element).html("Hide advanced settings");
	}
}

function toggleSugars(element) {
	var span = $(element);
	var sugars = $("#sugarTableHidden");
	
	if (sugars.is(":visible")) {
		span.children("#showText").css('display', 'block');
		span.children("#hideText").hide();
		sugars.hide();
	} else {
		span.children("#showText").hide();
		span.children("#hideText").css('display', 'block');
		sugars.show();
	}
}

function toggleSequence(element, event) {
	var link = $(element);
	var parent = $(element).parent();
	var sequence = parent.children('.sequence');
	
	if (sequence.is(":visible")) {
		link.hide();
		parent.children('.showLink').show();
		sequence.hide();
	} else {
		link.hide();
		parent.children('.hideLink').show();
		sequence.show();
	}

	safePreventEvent(event);     
    return false;
}

function safePreventEvent(e) {
	if (e.preventDefault) { 
		e.preventDefault()
	} else { 
		e.stop()
	};

    e.returnValue = false;
    e.stopPropagation();
}

function showDetails(element, event) {
	var next = $(element).parent().parent().parent().next('tbody.details');
	if (next.is(":visible")) {
		$(element).html("+");
		next.hide();
	} else {
		$(element).html("&ndash;");
		next.show();
	}

	safePreventEvent(event);
	return false;
}

var xmlhttp = new getXMLObject(); // xmlhttp holds the ajax object
function getXMLObject() {
	var xmlhttp;
	if (window.XMLHttpRequest) {// code for IE7+, Firefox, Chrome, Opera, Safari
		xmlhttp = new XMLHttpRequest();
	} else {// code for IE6, IE5
		xmlhttp = new ActiveXObject("Microsoft.XMLHTTP");
	}
	return xmlhttp; // Mandatory Statement returning the ajax object created
}

function initPage() {
	browserCheck();
	ajaxSessionRegister();
	$("select, input[type=radio], input[type=checkbox]").uniform();
	
	loadPageFromHash();
	
	$(window).hashchange(function () {
		if (this.location.hash != lastHash) {
			loadPageFromHash();
		}
	});
}

function browserCheck() {
	if (BrowserDetect.browser == "Explorer" || BrowserDetect.browser == "Mozilla"
		|| BrowserDetect.browser == "Netscape" || BrowserDetect.browser == "Konqueror") {
		$("#cboxClose").text('x');
		$.colorbox({ inline: true, href: "#browserError" });
	}
}

function loadPageFromHash() {
	var page = this.location.hash.replace('#!/', '');
	if (page == '') {
		page = 'prism';
	}
	var divs = ["help", "about", "json", "prism"];
	var works = false;
	for (var i = 0; i < divs.length; i++) {
		if (page == divs[i]) {
			works = true;
			break;
		}
	}
	if (!works) {
		return false;
	//	page = 'splash'; 
	} else {
		loadPage(page);		
	}
}

function loadPage(page) {
	// hide all divs
	var divs = ["helpDiv", "aboutDiv", "jsonDiv", "prismDiv"];
	for (var i = 0; i < divs.length; i++) {
		document.getElementById(divs[i]).style.display = "none";
	}
	
	useHashChangeEvent = false;
	var hash = '#!/';
	if (page == 'splash') {
		if (window.location.hash != '') {
			lastHash = hash;
			window.location.hash = hash;
		}
	} else {
		hash += page;
		lastHash = hash;
		window.location.hash = hash;
	}
	
	// smooth transitions with jquery
	// get an array of all divs, get the one showing, fade it out
	var shown;
	for (var i = 0; i < divs.length; i++) {
		if (document.getElementById(divs[i]).style.display != "block") shown = divs[i]; 
	}
	if (shown != null)
		$("#"+shown).fadeOut(700);
	$("#"+page+"Div").fadeIn(700);
}

function fileChange(id) {
	// simulate the default input[type=file] ui by getting file name
	var path = document.getElementById(id).value;
	if (path) {
		var startIndex = (path.indexOf('\\') >= 0 ? path.lastIndexOf('\\') : path.lastIndexOf('/'));
		var filename = path.substring(startIndex);
		if (filename.indexOf('\\') === 0 || filename.indexOf('/') === 0) {
			filename = filename.substring(1);
		}
		$("#"+id).parent().parent().parent().children(".inputFileName").html(filename);
	}
}

function ajaxSessionRegister() {
	var date = new Date();
	var month = date.getMonth() + 1;
	if (month.toString().length == 1) {
		month = '0' + month;
	}
	var day = date.getDate();
	if (day.toString().length == 1) {
		day = '0' + day;
	}
	var hour = date.getHours();
	if (hour.toString().length == 1) {
		hour = '0' + hour;
	}
	var mins = date.getMinutes();
	if (mins.toString().length == 1) {
		mins = '0' + mins;
	}
	
	var sessionID = date.getFullYear() + "" + month + "" + day + "-" + hour + "" + mins + "-"
		+ Math.floor(Math.random() * 1000000000);
	setSessionIDs(sessionID);
	
	params = "register=on&sessionID=" + sessionID;
	xmlhttp.open("POST", "PrismSessionSubmit", false);

	// Send the proper header information along with the request
	xmlhttp.setRequestHeader("Content-type", "application/x-www-form-urlencoded");
		//	"application/octet-stream");
			
	xmlhttp.onreadystatechange = sessionResponseHandler;
	xmlhttp.send(params);
}

var asks = 0;
function sessionResponseHandler() {
	if (xmlhttp.readyState == 4) {
		if (xmlhttp.status != 200) {
			alert("Session register error: " + xmlhttp.status + " : " + xmlhttp.statusText);
		} else {
			var response = xmlhttp.responseText;
			if (response.indexOf("response") == -1 && asks < 10) {
				window.setTimeout("ajaxSessionRegister();", 1000);
				asks++;
			}
		}
	}
}

function setSessionIDs(sessionID) {
	$("input[name=sessionID]").val(sessionID);
	$("input[name=sessionID]").html(sessionID);
	$("#sessionIDDiv").val(sessionID);
	$("#sessionIDDiv").html(sessionID);//convenience for debugging
}

function ajaxAcquireTaskInfo() {
	// generate timestamp
	date = new Date();
	stamp = date.getTime();
	var stampvalue = encodeURIComponent(stamp);
	
	// read session id
	var sessionvalue = $("#sessionIDDiv").val();
	
	// send
	xmlhttp.onreadystatechange = taskresponsehandler;
	xmlhttp.open("GET", "PrismSubmit" + "?sessionID=" + sessionvalue
			+ "&timestamp=" + stampvalue, true);
	xmlhttp.send();
}

var statusZero = 0;
var status500 = 0;
function taskresponsehandler() {
	if (xmlhttp.readyState == 4) {
		if (xmlhttp.status == 200) {
			try {
				document.getElementById("responseDiv").innerHTML = xmlhttp.responseText;
				var response = xmlhttp.responseText;
				if (response.indexOf("Job Done") == -1) {
					// Sleep then query again
					window.setTimeout("ajaxAcquireTaskInfo();", 1000);
				} else {
					// init sortable tables after task is done
					//	window.setTimeout("sorttable.init();", 100);							
				}
			} catch (err) {
				console.log(err);
				// Sleep then query again
				window.setTimeout("ajaxAcquireTaskInfo();", 1000);
			}

		} else if (xmlhttp.status == 0 && statusZero < 60) {
			window.setTimeout("ajaxAcquireTaskInfo();", 1000);
			statusZero++;
		} else if (xmlhttp.status == 500 && status500 < 60) {
			window.setTimeout("ajaxAcquireTaskInfo();", 1000);
			status500++;
		} else {
			if (xmlhttp.status == 0) {
				alert("Error: lost connection to the server for more than 60 seconds.");
			} else {
				alert(xmlhttp.status + " : " + xmlhttp.statusText);
			}
		}
	}
}

function submitCheck() {
	var filename1 = $("#genomeFile").val();
	if (filename1 == null || filename1 == "") {
		alert("Genome: Please select a FASTA or GenBank file!");
		return false;
	}
	{
		window.setTimeout("disableForm('genomeForm');", 100);
		var session = $("#sessionIDDiv").val();
		$("#buttonHolder").html("<p class='thanks'>Your results will be available at <br><a href='/prism/tasks/" 
				+ session + "/'>/prism/tasks/" + session + "/</a>.	</p>").show();
		window.setTimeout("ajaxAcquireTaskInfo();", 1000);
		return true;
	}
}

function jsonSubmitCheck() {
	var filename = $("#jsonFile").val();
	if (filename == null || filename == "") {
		alert("Load saved results: Please select a JSON file!");
		return false;
	}
	{
		window.setTimeout("disableForm('jsonForm');", 100);
		var session = $("#sessionIDDiv").val();
		$("#jsonButtonHolder").html("<p class='thanks'>Your results will be available at <br><a href='/prism/tasks/" 
				+ session + "/'>/prism/tasks/" + session + "/</a>.	</p>").show();
		window.setTimeout("ajaxAcquireJsonTaskInfo();", 1000);
		return true;
	}
}

function ajaxAcquireJsonTaskInfo() {
	// generate timestamp
	date = new Date();
	stamp = date.getTime();
	var stampvalue = encodeURIComponent(stamp);
	
	// read session id
	var sessionvalue = $("#sessionIDDiv").val();
	
	// send
	xmlhttp.onreadystatechange = jsonresponsehandler;
	xmlhttp.open("GET", "PrismJsonSubmit" + "?sessionID=" + sessionvalue
			+ "&timestamp=" + stampvalue, true);
	xmlhttp.send();
}

var jsonStatusZero = 0;
var jsonStatus500 = 0;
function jsonresponsehandler() {
	if (xmlhttp.readyState == 4) {
		if (xmlhttp.status == 200) {
			try {
				document.getElementById("jsonResponseDiv").innerHTML = xmlhttp.responseText;
				var response = xmlhttp.responseText;
				if (response.indexOf("Job Done") == -1) {
					// Sleep then query again
					window.setTimeout("ajaxAcquireJsonTaskInfo();", 1000);
				}
			} catch (err) {
				console.log(err);
				// Sleep then query again
				window.setTimeout("ajaxAcquireJsonTaskInfo();", 1000);
			}

		} else if (xmlhttp.status == 0 && jsonStatusZero < 60) {
			window.setTimeout("ajaxAcquireJsonTaskInfo();", 1000);
			jsonStatusZero++;
		} else if (xmlhttp.status == 500 && jsonStatus500 < 60) {
			window.setTimeout("ajaxAcquireJsonTaskInfo();", 1000);
			jsonStatus500++;
		} else {
			if (xmlhttp.status == 0) {
				alert("Error: lost connection to the server for more than 60 seconds.");
			} else {
				alert(xmlhttp.status + " : " + xmlhttp.statusText);
			}
		}
	}
}

function disableForm(id) {
	var form = $("#"+id);
	var $inputs = $("#" + id + " :input");
	for (var i = 0; i < $inputs.length; i++) {
		var tempobj = $inputs[i];
		if (tempobj.type != null
				&& (tempobj.type.toLowerCase() == "submit" || tempobj.type
						.toLowerCase() == "reset")) {
			tempobj.parentElement.style.display = "none";
		} else {
			tempobj.disabled = true;
		}
	}
}

function numbersOnly(myfield, e) {
	var key;
	var keychar;
	if (window.event)
		key = window.event.keyCode;
	else if (e)
		key = e.which;
	else
		return true;
	keychar = String.fromCharCode(key);

	// numbers
	if ((("0123456789").indexOf(keychar) > -1))
		return true;
	// decimal point
	else if (keychar == '.' && myfield.value.indexOf('.') < 0)
		return true;
	else
		return false;
}

function integersOnly(myfield, e) {
	var key;
	var keychar;
	if (window.event)
		key = window.event.keyCode;
	else if (e)
		key = e.which;
	else
		return true;
	keychar = String.fromCharCode(key);

	// numbers
	if ((("0123456789").indexOf(keychar) > -1))
		return true;
	else
		return false;
}

function endsWith(str, suffix) {
	return str.length >= suffix.length
			&& str.substr(str.length - suffix.length) == suffix;
}

var BrowserDetect = {
		init: function () {
			this.browser = this.searchString(this.dataBrowser) || "An unknown browser";
			this.version = this.searchVersion(navigator.userAgent)
				|| this.searchVersion(navigator.appVersion)
				|| "an unknown version";
			this.OS = this.searchString(this.dataOS) || "an unknown OS";
		},
		searchString: function (data) {
			for (var i=0;i<data.length;i++)	{
				var dataString = data[i].string;
				var dataProp = data[i].prop;
				this.versionSearchString = data[i].versionSearch || data[i].identity;
				if (dataString) {
					if (dataString.indexOf(data[i].subString) != -1)
						return data[i].identity;
				}
				else if (dataProp)
					return data[i].identity;
			}
		},
		searchVersion: function (dataString) {
			var index = dataString.indexOf(this.versionSearchString);
			if (index == -1) return;
			return parseFloat(dataString.substring(index+this.versionSearchString.length+1));
		},
		dataBrowser: [
			{
				string: navigator.userAgent,
				subString: "Chrome",
				identity: "Chrome"
			},
			{ 	string: navigator.userAgent,
				subString: "OmniWeb",
				versionSearch: "OmniWeb/",
				identity: "OmniWeb"
			},
			{
				string: navigator.vendor,
				subString: "Apple",
				identity: "Safari",
				versionSearch: "Version"
			},
			{
				prop: window.opera,
				identity: "Opera",
				versionSearch: "Version"
			},
			{
				string: navigator.vendor,
				subString: "iCab",
				identity: "iCab"
			},
			{
				string: navigator.vendor,
				subString: "KDE",
				identity: "Konqueror"
			},
			{
				string: navigator.userAgent,
				subString: "Firefox",
				identity: "Firefox"
			},
			{
				string: navigator.vendor,
				subString: "Camino",
				identity: "Camino"
			},
			{		// for newer Netscapes (6+)
				string: navigator.userAgent,
				subString: "Netscape",
				identity: "Netscape"
			},
			{
				string: navigator.userAgent,
				subString: "MSIE",
				identity: "Explorer",
				versionSearch: "MSIE"
			},
			{
				string: navigator.userAgent,
				subString: "Gecko",
				identity: "Mozilla",
				versionSearch: "rv"
			},
			{ 		// for older Netscapes (4-)
				string: navigator.userAgent,
				subString: "Mozilla",
				identity: "Netscape",
				versionSearch: "Mozilla"
			}
		],
		dataOS : [
			{
				string: navigator.platform,
				subString: "Win",
				identity: "Windows"
			},
			{
				string: navigator.platform,
				subString: "Mac",
				identity: "Mac"
			},
			{
				   string: navigator.userAgent,
				   subString: "iPhone",
				   identity: "iPhone/iPod"
		    },
			{
				string: navigator.platform,
				subString: "Linux",
				identity: "Linux"
			}
		]

	};
	BrowserDetect.init();