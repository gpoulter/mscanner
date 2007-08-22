function readFile(fname) {
   if (window.ActiveXObject) {
      return readFileIE(fname);
   } else if (netscape) {
      return readFileNetscape(fname);
   } else {
      alert("Can only read files from IE or Netscape/Mozilla/Firefox");
   }
}

function writeFile(fname, text) {
   if (window.ActiveXObject) {
      writeFileIE(fname, text);
   } else if (netscape) {
      writeFileNetscape(fname, text);
   } else {
      alert("Can only save files from IE or Netscape/Mozilla/Firefox");
   }
}

function readFileIE(fname, text) {
   fso  = new window.ActiveXObject("Scripting.FileSystemObject");
   os = fso.OpenTextFile(fname, 1);
   result = os.ReadAll();
   os.Close();
   return result;
}

function writeFileIE(fname, text) {
   fso  = new window.ActiveXObject("Scripting.FileSystemObject");
   os = fso.OpenTextFile(fname, 2, true);
   os.Write(text);
   os.Close();
}

function readFileNetscape(fname) {
   try {
      netscape.security.PrivilegeManager.enablePrivilege("UniversalXPConnect");
   } catch (e) {
      alert("Permission to read file was denied.");
   }
   var file = Components.classes["@mozilla.org/file/local;1"].createInstance(Components.interfaces.nsILocalFile);
   file.initWithPath(fname);
   if (file.exists() == false) {
      alert("File does not exist");
   }
   var is = Components.classes["@mozilla.org/network/file-input-stream;1"].createInstance(Components.interfaces.nsIFileInputStream );
   is.init( file,0x01, 00004, null);
   var sis = Components.classes["@mozilla.org/scriptableinputstream;1"].createInstance(Components.interfaces.nsIScriptableInputStream );
   sis.init(is);
   return sis.read(sis.available());
}

function writeFileNetscape(fname, text) {
   try {
      netscape.security.PrivilegeManager.enablePrivilege("UniversalXPConnect");
   } catch (e) {
      alert("Permission to save file was denied.");
   }
   var file = Components.classes["@mozilla.org/file/local;1"].createInstance(Components.interfaces.nsILocalFile);
   file.initWithPath(fname);
   if (file.exists() == false) {
      file.create(Components.interfaces.nsIFile.NORMAL_FILE_TYPE, 420);
   } 
   var outputStream = Components.classes["@mozilla.org/network/file-output-stream;1"].createInstance(Components.interfaces.nsIFileOutputStream);
   outputStream.init(file, 0x04 | 0x08 | 0x20, 420, 0);
   var result = outputStream.write(text, text.length);
   outputStream.close();
}

/* Load XML document, calls back when ready. */
// http://www.quirksmode.org/dom/importxml.html 
function importXML(fname, callback)
{
   if (document.implementation && document.implementation.createDocument)
   {
      xmlDoc = document.implementation.createDocument("", "", null);
      xmlDoc.onload = callback;
   }
   else if (window.ActiveXObject)
   {
      xmlDoc = new ActiveXObject("Microsoft.XMLDOM");
      xmlDoc.onreadystatechange = function () {
	 if (xmlDoc.readyState == 4) callback(xmlDoc);
      };
   }
   else
   {
      alert('Your browser can\'t handle this script');
      return;
   }
   xmlDoc.load(filename);
}
