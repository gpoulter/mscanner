/* Store all global variables in global object "g" */
var g = new Object();

/* Things to do as soon as the HTML is loaded */
window.onload = function() {
   g.outdir = "output";
   g.batchid = "";
   // Set up jsolait XML-RPC
   //jsolait.baseURI = "/htdocs/script/jsolait";
   g.xmlrpc = imprt("xmlrpc");
   //var service = "http://jonbn.med.uct.ac.za/cgi-bin/mscanner.py"
   var service = "http://medline.stanford.edu/cgi-bin/mscanner.py"
   var methods = ["getStatus","listBatches","deleteBatch","query","validate"]
   g.service = new g.xmlrpc.ServiceProxy(service, methods);
   // Set up sliding help effects
   g.fx = new Object();
   var names = [
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
   getStatus();
}

/* Hide an element */
function hideBlock(block) {
   $(block).style.display = "none";
}

/* Unhide an element */
function showBlock(block) {
   $(block).style.display = "block";
}

/* Delete batch output */
function deleteBatch(batchid) {
   var status_callback = function(result, err) {
      if (err != null) {
         alert("XMLRPC Error in getStatus():\n\n" + err.faultString);
         return
      }
      if (result.status == "busy") {
         alert("Cannot delete since busy processing another batch");
         return;
      }
      var delete_callback = function(result, err) {
         if (err != null) {
            alert("XMLRPC Error in deleteBatch():\n\n" + err.faultString);
            return
         }
         $("status").innerHTML = ""
         graft($("status"),['h2.heading', "Batch Deleted"])
         graft($("status"),[
            'div.body',
            ['p.para',
             "Batch ID " + batchid + " has been deleted."
            ]
         ])
         $("status").innerHTML = $("status").innerHTML
         showBlock("status");
         listBatches();
      }
      g.service.deleteBatch(batchid, delete_callback);
   }
   g.service.getStatus("", status_callback);
}

/* Start a query, returning batch id */
function query() {
   if (g.rpcbusy) {
      alert("A request has already been sent - " +
            "just waiting for reply. If you think the " +
            "reply should be here by now, refresh the page.");
      return;
   }
   var query_callback = function(result, err) {
      g.rpcbusy = false;
      if (err != null) {
         alert("Cannot process query request: \n\n" + err.faultString);
         return;
      }
      g.batchid = result;
      setTimeout("getStatus()", 1000);
   }
   var positives = parsePMIDs($("positives").value);
   var pseudocount = parseFloat($("pseudocount").value);
   var limit = parseInt($("limit").value);
   var threshold = parseFloat($("threshold").value);
   g.rpcbusy = true;
   g.service.query(positives, pseudocount, limit, threshold, query_callback);
}

/* Start validation, returning batch id */
function validate() {
   if (g.rpcbusy) {
      alert("A request has already been sent - " +
            "just waiting for reply. If you think the " +
            "reply should be here by now, refresh the page.");
      return;
   }
   var valid_callback = function(result, err) {
      g.rpcbusy = false;
      if (err != null) {
         alert("Cannot process validation request: \n\n" + err.faultString);
         return;
      }
      g.batchid = result;
      setTimeout("getStatus()", 1000);
   }
   var positives = parsePMIDs($("positives").value);
   var negatives = parseInt($("negatives").value);
   var nfold = parseInt($("nfold").value);
   var pseudocount = parseFloat($("pseudocount").value);
   g.rpcbusy = true;
   g.service.validate(positives, negatives, nfold, pseudocount, valid_callback);
}

/* Print whether scanner is busy or done with this file.  If busy,
call itself again in 1 second.  */
function getStatus() {
   var status_callback = function(result, err) {
      if (err != null) {
         alert("XMLRPC Error in getStatus(): \n\n" + err.faultString);
         return;
      }
      if (g.batchid == "") {
         g.batchid = result.batchid;
      }
      if (result.status == "busy") {
         $("status").innerHTML = ""
         graft($("status"),['h2.heading', "Busy Scanning!"])
         graft($("status"),[
            'div.body',
            ['p.para',
             "Batch ID " + g.batchid + " is busy running."
            ],
            ['p.para',
             "Note that queries may take up to and hour, and validation takes "+
             "time proportional to the number of articles involved. "+
             "However, the output will be placed in ",
             ['a',
              {href: g.outdir + "/" + g.batchid + "/index.html"},
              g.batchid
             ],
             " which you can bookmark and come back to later for the results."
            ],
            ['p.para',
             "It has been " + Math.floor(result.elapsed/60) + " minutes and " + 
             Math.floor(result.elapsed % 60) + " seconds since scanning started, and "+
             result.progress + " steps out of " + result.total + " have been completed."
            ]
         ])
         $("status").innerHTML = $("status").innerHTML
         showBlock("status");
         setTimeout("getStatus()", 4000);
      }
      else if (result.status == "done") {
         $("status").innerHTML = ""
         graft($("status"),['h2.heading', "Scanning Complete!"])
         graft($("status"),[
            'div.body',
            ['p.para',
             "Results are available for batch ID ",
             ['a',
              {href: g.outdir + "/" + g.batchid + "/index.html"},
              ['b',
               g.batchid
              ]
             ]
            ]
         ])
         $("status").innerHTML = $("status").innerHTML
         showBlock("status");
         listBatches();
      }
      else if (result.status == "notfound") {
         $("status").innerHTML = ""
         alert("Batch " + g.batchid + " has failed!\n" +
               "E-mail xxxxxxxxxxxxxxxxxxxxxxxx with questions.\n" +
               "Script output is:\n" + result.lastlog + "\n"
              );
      }
   }
   g.service.getStatus(g.batchid, status_callback);
}

/* List available results */
function listBatches() {
   var listing_callback = function(result, err) {
      if (err != null) {
         alert("XMLRPC Error in listBatches(): \n\n" + err.faultString);
         return;
      }
      $("batches").innerHTML = ""
      if (result.length == 0) {
         hideBlock("batches");
      } else {
         $("batches").innerHTML = "";
         graft($("batches"), ['h2.heading', "Results" ])
         for (var x = 0; x < result.length; x++) {
            var batch = result[x];
            graft($("batches"),[
               'p.item',
               ['button.delete',
                {onclick: "deleteBatch('" + batch +"');"},
                "Delete"
               ],
               ['a',
                {href: g.outdir + "/" + batch + "/index.html"},
                ['b',
                 batchIdToDate(batch)
                ]
               ]
            ])
         }
         $("batches").innerHTML = $("batches").innerHTML
         showBlock("batches");
      }
   }
   g.service.listBatches(listing_callback);
}

/* Convert 'YYYYmmdd-HHMMSS' to 'YYYY/mm/dd HH:MM:SS' */
function batchIdToDate(batchid) {
   var b = batchid;
   if (b.match(/^\d{8}-\d{6}$/)) {
      return "".concat(b.substring(0,4),"/",b.substring(4,6),"/",b.substring(6,8)," ",
                       b.substring(9,11),":",b.substring(11,13),":",b.substring(13,15));
   } else {
      return batchid;
   }
}

/* Convert a string of newline-separated PubMed IDs into a list of integers */
function parsePMIDs(pmidstr) {
   try {
      spmids = pmidstr.split("\n");
      ipmids = [];
      for(var i = 0; i < spmids.length; i++) {
         ipmids.push(parseInt(spmids[i]));
      }
      return ipmids;
   }
   catch(e) {
      alert("Bad formatting in PMID list: " + e);
      return null;
   }
}