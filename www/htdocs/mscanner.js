/* Store all global variables in global object "g" */
var g = new Object();

/* Things to do as soon as the HTML is loaded */
window.onload = function() {
   g.outdir = "output";
   g.batchid = "";
   // Set up jsolait XML-RPC
   //jsolait.baseURI = "/htdocs/script/jsolait";
   g.xmlrpc = imprt("xmlrpc");
   g.service = new g.xmlrpc.ServiceProxy(
      "http://jonbn.med.uct.ac.za/cgi-bin/mscanner.py",
      //"http://medline.stanford.edu/cgi-bin/mscanner.py",
      ["getStatus","listBatches","deleteBatch","query","validate"]);
   // Set up sliding help effects
   g.fx = new Object();
   names = [
      "positives_help",
      "pseudocount_help",
      "limit_help",
      "threshold_help",
      "negatives_help",
      "nfold_help"
   ];
   for (var i = 0; i < names.length; i++) {
      g.fx[names[i]] = new fx.Combo( names[i] , {height: true, opacity: true, duration: 300} );
      g.fx[names[i]].hide();
   }
   listBatches();
}

/* SXOOP templating of destination.innerHTML using source.value */
function doTemplate(destination, source, params) {
   $(destination).innerHTML = SXOOP.template.parse($(source).value, params);
}

/* Hide an element and all children */
function hideBlock(block) {
   $(block).style.display = "none";
}

/* Unhide an element and all children *
function showBlock(block) {
   $(block).style.display = "block";
}

/* Delete batch output */
function deleteBatch(batchid) {
   result = g.service.getStatus("");
   if (result.status == "busy") {
      alert("Cannot delete: busy processing another batch");
      return;
   }
   try {
      g.service.deleteBatch(batchid);
   } catch(e) {
      alert("Error: \n\n" + e.message);      
      return;
   }
   doTemplate("status", "t_deleted", batchid);
}

/* Start a query, returning batch id */
function query() {
   //alert($("positives"));
   positives = $("positives").value;
   pseudocount = parseFloat($("pseudocount").value);
   limit = parseInt($("limit").value);
   threshold = parseFloat($("threshold").value);
   try {
      batchid = g.service.query(positives, pseudocount, limit, threshold);
   } catch(e) {
      alert("Error: \n\n" + e.message);
      return;
   }
   g.batchid = batchid
   setTimeout("getStatus()", 1000);
}

/* Start validation, returning batch id */
function validate() {
   positives = $("positives").value;
   negatives = parseInt($("negatives").value);
   nfold = parseInt($("nfold").value);
   pseudocount = parseFloat($("pseudocount").value);
   try {
      batchid = g.service.validate(positives, negatives, nfold, pseudocount);
   } catch(e) {
      alert("Error: \n\n" + e.message);
      return;
   }
   g.batchid = batchid
   setTimeout("getStatus()", 1000);
}

/* Print whether scanner is busy or done with this file.  If busy,
call itself again in 1 second.  */
function getStatus() {
   try {
      result = g.service.getStatus(g.batchid);
      if (g.batchid == "") {
         g.batchid = result.batchid;
      }
   } catch(e) {
      alert("Error: \n\n" + e.message);
      return;
   }
   if (result.status == "busy") {
      doTemplate("status", "t_busy", g.batchid);
      g.started = result.started * 1000;
      setTimeout("getStatus()", 1000);
   }
   else if (result.status == "done") {
      doTemplate("status", "t_done", g.batchid);
      listBatches();
   }
   else if (result.status == "notfound") {
      alert("Batch " + g.batchid + " has failed!\n" +
            "E-mail xxxxxxxxxxxxxxxxxxxxxxxx with questions.\n" +
            "Script output is:\n" + result.lastlog + "\n"
             );
   }
}

/* List available results */
function listBatches() {
   try {
      batches = g.service.listBatches();
   } catch(e) {
      alert("Error: \n\n" + e.message);
      return;
   }
   doTemplate("batches", "t_batches", batches);
}

/* Convert YYYYmmdd-HHMMSS -> YYYY/mm/dd HH:MM:SS */
function batchIdToDate(batchid) {
   b = batchid;
   return "".concat(b.substr(0,4),"/",b.substr(4,2),"/",b.substr(6,2)," ",
                    b.substr(9,2),":",b.substr(11,2),":",b.substr(13,2));
}