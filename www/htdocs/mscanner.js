/* Store all global variables in global object "g" */
var g = new Object();

/* Things to do as soon as the HTML is loaded */
window.onload = function() {
   g.outdir = "output";
   g.batchid = "";
   // Set up jsolait XML-RPC
   //jsolait.baseURI = "/htdocs/script/jsolait";
   g.xmlrpc = imprt("xmlrpc");
   var service = "http://jonbn.med.uct.ac.za/cgi-bin/mscanner.py"
   //var service = "http://mscanner.stanford.edu/cgi-bin/mscanner.py"
   //var service = "http://medline.stanford.edu/cgi-bin/mscanner.py"
   var methods = ["getStatus","listBatches","deleteBatch","query","validate"]
   g.service = new g.xmlrpc.ServiceProxy(service, methods);
   // Set up sliding help effects
   g.fx = new Object();
   var names = [
      "title_help",
      "code_help",
      "positives_help",
      "pseudocount_help",
      "limit_help",
      "threshold_help",
      "nfold_help",
      "negatives_help",
      "alpha_help"
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
function deleteBatch(batchid, delcode) {
   var status_callback = function(result, err) {
      if (err != null) {
         alert("XMLRPC Error in getStatus():\n\n" + err.faultString);
         return
      }
      if (result.status == "busy") {
         alert("Cannot delete since busy processing another batch");
         return;
      }
      /* If status is ok for deletion, do it! */
      var delete_callback = function(result, err) {
         if (err != null) {
            alert("XMLRPC Error in deleteBatch():\n\n" + err.faultString);
         } else {
            alert("Results for Batch ID " + batchid + " have been deleted!");
         }
         listBatches();
      }
      g.service.deleteBatch(batchid, delcode, delete_callback);
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
         g.batchid = null;
         return;
      }
      setTimeout("getStatus()", 1000);
   }
   var batchid = $("title")
   var code = $("code")
   var positives = parsePMIDs($("positives").value);
   var pseudocount = parseFloat($("pseudocount").value);
   var limit = parseInt($("limit").value);
   var threshold = parseFloat($("threshold").value);
   g.rpcbusy = true;
   g.batchid = batchid;
   g.service.query(batchid, code, positives, pseudocount, limit, threshold, query_callback);
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
         g.batchid = null;
         return;
      }
      setTimeout("getStatus()", 1000);
   }
   var batchid = $("title")
   var code = $("code")
   var positives = parsePMIDs($("positives").value);
   var negatives = parseInt($("negatives").value);
   var nfold = parseInt($("nfold").value);
   var pseudocount = parseFloat($("pseudocount").value);
   var alpha = parseFloat($("alpha"));
   g.rpcbusy = true;
   g.batchid = batchid;
   g.service.validate(batchid, code, positives, negatives, nfold, pseudocount, alpha, valid_callback);
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
             "Batch ID " + g.batchid + " is busy running. It has been " +
             Math.floor(result.elapsed / 60) +
             " minutes and " +
             Math.floor(result.elapsed % 60) +
             " seconds since scanning started, and "+
             result.progress + " steps out of " + result.total + " have been completed."
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
             "If instead you wish to be notified when the scanner is available," +
             "leave your e-mail address (nobody will see it, and it will be " +
             "deleted after sending notification):",
             ['input',
              {class:"field", type:"text", size:"30", value:"", id:"email"},
             ],
             ['button',
              {onclick:"addMailer($('email'))"},
              "submit"
             ]
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
         graft($("batches"),
               ['h2.heading',
                "Results"
               ]
              )
         graft($("batches"),
               ['p.para',
                "To delete your batch, enter the deletion code you chose:",
                ['input',
                 {class:"field", type:"text", size:"8", value:"", id:"delete_code"}
                ],
               ]
              )
         for (var x = 0; x < result.length; x++) {
            var batch = result[x];
            graft($("batches"),[
               'p.item',
               ['button.delete',
                {onclick: "deleteBatch('" + batch +"', $('delete_code'));"},
                "Delete"
               ],
               ['a',
                {href: g.outdir + "/" + batch + "/index.html"},
                ['b',
                 batch
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

/* Add an e-mail address to the list to be alerted when the scanner is available */
function addMailer(email) {
   var mailer_callback = function(result, err) {
      if (err != null) {
         alert("Failed to add e-mail address: \n\n" + err.faultString);
         return;
      }
   }
   g.service.addMailer(email, query_callback);
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
