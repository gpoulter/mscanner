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
   var result = g.service.getStatus("");
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

/* Start a query, returning batch id */
function query() {
   //alert($("positives"));
   var positives = $("positives").value;
   var pseudocount = parseFloat($("pseudocount").value);
   var limit = parseInt($("limit").value);
   var threshold = parseFloat($("threshold").value);
   try {
      g.batchid = g.service.query(positives, pseudocount, limit, threshold);
      setTimeout("getStatus()", 1000);
   } catch(e) {
      alert("Error: \n\n" + e.message);
      return;
   }
}

/* Start validation, returning batch id */
function validate() {
   var positives = $("positives").value;
   var negatives = parseInt($("negatives").value);
   var nfold = parseInt($("nfold").value);
   var pseudocount = parseFloat($("pseudocount").value);
   try {
      g.batchid = g.service.validate(positives, negatives, nfold, pseudocount);
      setTimeout("getStatus()", 1000);
   } catch(e) {
      alert("Error: \n\n" + e.message);
      return;
   }
}

/* Print whether scanner is busy or done with this file.  If busy,
call itself again in 1 second.  */
function getStatus() {
   try {
      var result = g.service.getStatus(g.batchid);
      if (g.batchid == "") {
         g.batchid = result.batchid;
      }
   } catch(e) {
      alert("Error: \n\n" + e.message);
      return;
   }
   if (result.status == "busy") {
      $("status").innerHTML = ""
      graft($("status"),['h2.heading', "Busy Scanning!"])
      graft($("status"),[
         'div.body',
         ['p.para',
          "Batch ID " + g.batchid + " is still running."
         ],
         ['p.para',
          "Note that queries may take up to 20 minutes, and validation "+
          "time proportional to the number of articles being processed. "+
          "The output will be placed in ",
          ['a',
           {href: g.outdir + "/" + g.batchid + "/index.html"},
           g.batchid
          ],
          " so you can save the link and come back later if you please."
         ],
         ['p.para',
          "It has been " + result.elapsed + " seconds since scanning started, and "+
          result.progress + " steps out of " + result.total + " have been completed."
         ]
      ])
      $("status").innerHTML = $("status").innerHTML
      showBlock("status");
      setTimeout("getStatus()", 1000);
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
      alert("Batch " + g.batchid + " has failed!\n" +
            "E-mail xxxxxxxxxxxxxxxxxxxxxxxx with questions.\n" +
            "Script output is:\n" + result.lastlog + "\n"
             );
   }
}

/* List available results */
function listBatches() {
   try {
      var batches = g.service.listBatches();
   } catch(e) {
      alert("Error: \n\n" + e.message);
      return;
   }
   $("batches").innerHTML = ""
   if (batches.length == 0) {
      hideBlock("batches");
   } else {
      $("batches").innerHTML = "";
      graft($("batches"), ['h2.heading', "Results" ])
      for (var x = 0; x < batches.length; x++) {
         var batch = batches[x];
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

/* Convert YYYYmmdd-HHMMSS -> YYYY/mm/dd HH:MM:SS */
function batchIdToDate(batchid) {
   var b = batchid;
   return "".concat(b.substr(0,4),"/",b.substr(4,2),"/",b.substr(6,2)," ",
                    b.substr(9,2),":",b.substr(11,2),":",b.substr(13,2));
}

/* Convert a string of newline-separated PubMed IDs into a list of integers */
function parsePMIDs(pmidstr) {
   
}