/*
Cells in each row are:
0 - Classification
1 - Rank
2 - Score
3 - PMID
4 - Year
5 - Expand Author
6 - Expand Abstract
7 - Title
8 - ISSN
*/

/************************ HELPER FUNCTIONS ***********************/

/* Return whether event was a right click */
function is_rightclick(event) {
   if (!event) var event = window.event;
   if (event.which)
      return (event.which == 3);
   else if (event.button)
      return (event.button == 2);
}

/* Return target of the event */
function event_target(event) {
   if (!event) var event = window.event;
   var targ = null;
   if (event.target) 
      targ = event.target; // Netscape
   else if (event.srcElement)
      targ = event.srcElement; // Microsoft
   if (targ.nodeType == 3) // defeat Safari bug
      targ = targ.parentNode;
   return targ;
}  

/* return a bound method, can be used to create event handlers
that don't switch to the context of the HTML element */
function bind(obj, method) {
   return function() { return method.apply(obj, arguments); }
}

/* Use onclick="return noenter()" to prevent default form submission 
action when enter key is pressed in an input field  */
function noenter() {
  return !(window.event && window.event.keyCode == 13); 
}

/********************* OPENING RESULTS IN PUBMED ***************************/

/* Open positive documents in PubMed 
BUG: Internet Explorer barfs on long query strings!

all=true: Open all visible citations in PubMed
all=false: Open only manually-marked relevant citations in PubMed
*/
function openPubMed(all) {
   var rows = document.getElementById("citations").rows;
   var qstring = "http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&amp;db=pubmed&amp;list_uids=";
   var count = 0;
   for (var i = 0; i < rows.length; i++) {
      if (rows[i].className == "main") {
         var pmid = rows[i].id;
         if (getTag(pmid) == 2 ||
            (all == true && rows[i].style.display != "none"))
         {
            if (count > 0) {
               qstring += ",";
            }
            qstring += pmid;
            count++;
         }
      }
   }
   if(count == 0) {
      alert("Cannot open PubMed as there are no PubMed IDs to open");
   } else {
      window.open(qstring);
   }
}

/* Adds events for opening citations in PubMed */
function initPubMedOpening() {
   $('open_pubmed_relevant').onclick = function() {
      openPubMed(false); }
   $('open_pubmed_all').onclick = function() {
      openPubMed(true);  }
}

/******************* LOADING/SAVING OF CITATION TAGS ********************/

/*
0 = Unclassified (gray)
1 = Negative (blue)
2 = Positive (red)
*/

/* Read classifications from disk */
function loadTags(fname) 
{
   try {
      var text = readFile(fname)
   } catch (err) {
      alert("Failed to read classifications from file: " + err);
      return;
   }
   var lines = text.split("\n");
   for (var i = 0; i < lines.length; i++) {
      var values = lines[i].split(",");
      var pmid = values[0];
      var score = values[1];
      var classification = values[2];
      switch(classification) {
         case " ": setTag(pmid, 0); break;
         case "0": setTag(pmid, 1); break;
         case "1": setTag(pmid, 2); break;
      }
   }
   alert("Successfully loaded classifications");
}

/* Save classifications to disk */
function saveTags(fname) {
   var text = "";
   var rows = document.getElementById("citations").rows;
   for (var i = 0; i < rows.length; i++) {
      if (rows[i].className == "main") {
         var pmid = rows[i].id;
         var score = rows[i].cells[2].innerHTML;
         var line = pmid+","+score+","
         switch(getTag(pmid)) {
            case 0: line += " \n"; break;
            case 1: line += "0\n"; break;
            case 2: line += "1\n"; break;
         }
         text += line;
      }
   }
   try {
      writeFile(fname, text);
      alert("Successfully saved classifications");
   } catch (err) {
      alert("Failed to save classifications: " + err);
   }
}

/* Set up the saving/loading of classifications */
function initTagSaving() {
   if (areWeLocal() == false) {
      // hide ability to save/load citations for non-local pages
      $("save_load_div").style.display = "none";
   } else {
      // else hide warning and hook up the events
      $("no_save_load_div").style.display = "none";
      $('save_tags').onclick = function() {
         saveTags(fileURLasPath($('save_target').href));
      }
      $('load_tags').onclick = function() {
         loadTags(fileURLasPath($('save_target').href));
      }
   }
}

/******************** TAGGING OF CITATIONS ************************/

/* Get classification tag (0,1,2) */
function getTag(pmid) {
   var color = $(pmid).cells[0].style.backgroundColor
   switch (color) {
      case "": return 0;
      case "blue": return 1;
      case "red": return 2;
      default: throw "Invalid tag color " + color;
   } 
}

/* Set classification tag (value may be 0, 1 or 2) */
function setTag(pmid, value) {
   var row = $(pmid)
   var tagstyle = row.cells[0].style
   var rowstyle = row.style
   switch (value) {
   case 0:
      tagstyle.backgroundColor = "";
      rowstyle.backgroundColor = "";
      break;
   case 1:
      tagstyle.backgroundColor = "blue"; 
      rowstyle.backgroundColor = "#DDDDFF";
      break;
   case 2:
      tagstyle.backgroundColor = "red";
      rowstyle.backgroundColor = "#FFDDDD";
      break;
   default: throw "Invalid tag value " + value;
   } 
}


/******************* CITATION TABLE EVENTS **************************/

/* Hide all expanded author/abstract rows in the citations table */
function hideTableRows() {
   var rows = document.getElementById("citations").rows;
   // row 0 is the heading
   for (i = 1; i < rows.length; i+=3) {
      rows[i].style.display = "";
      rows[i+1].style.display = "none";            
      rows[i+2].style.display = "none";
   }
}

/* Adds the events for classification/author/abstract display */
function initTableEvents() {

   hideTableRows()

   /* Cycles cell forward on left clicks, backward on right */
   function onclick_classify(e) {
      if (!e) var e = window.event;
      // Target is cell, with parent row who has PMID as its id
      var pmid = this.parentNode.id;
      if (is_rightclick(e)) {
         setTag(pmid, (getTag(pmid)+1) % 3);
      } else {
         setTag(pmid, (getTag(pmid)+2) % 3);
      }
      return false; // suppress default click action
   }
   
   /* Author row is sibling of target's parent row. */
   function onclick_author() {
      toggle(this.parentNode.nextSibling);
   }

   /* Abstract row is second sibling of target's parent row. */
   function onclick_abstract() {
      toggle(this.parentNode.nextSibling.nextSibling);
   }
   
   /* Loop over rows to add event handlers */
   var rows = document.getElementById("citations").rows;
   // row 0 is the heading
   for (var i = 1; i < rows.length; i+=3) {
      rows[i].cells[0].onclick = onclick_classify;
      rows[i].cells[5].onclick = onclick_author;
      rows[i].cells[6].onclick = onclick_abstract;
   }
}


/****************** APPENDING TO THE CITATION TABLE ************************/

/* 
Read citations from disk and append to the table
If append is false, we prepend the data to the table instead.
*/
function loadCitations(fname) {
   var text = readFile(fname);
   var table = text.match(/<[t]body>[^\v]*?<\/tbody>/);
   var blob = table[0].substr(7,table[0].length-15);
   var tbody = $("citations").getElementsByTagName("tbody")[0]
   tbody.innerHTML += blob;
   /* new HTML so no events are linked */
   initTableEvents();
}

/* Note: We can only append result files if innerHTML can be set on tables
(works at least in mozilla browsers) and we are a local file. */

appendable_files = [] // list of files to append
function initResultAppending() {
   /* hide the feature if it is not supported. toSource detects Netscape. */
   if (areWeLocal() == false || (!"a".toSource)) {
      $("append_results_div").style.display = "none";
      return;
   }
   /* make a list of files to load */
   nodes = $("appendable_files").getElementsByTagName("a");
   for(var i = 0; i < nodes.length; i++) {
      appendable_files.push(nodes[i].href);
   }
   function setupNextAppend() {
      if (appendable_files.length == 0) {
         $("next_to_append").href = "#";
         $("next_to_append").innerHTML = "(no more files to append)";
      } else {
         function fn(url) { return url.substr(url.lastIndexOf("/")+1); }
         $("next_to_append").href = appendable_files.shift();
         $("next_to_append").innerHTML = fn($("next_to_append").href);
      }
   }
   $('append_results').onclick = function() {
      if ($("next_to_append").href != "#") {
         loadCitations(fileURLasPath($("next_to_append").href));
         setupNextAppend();
      }
   }
   setupNextAppend();
}

/*********************** FILTERING OF CITATIONS ************************/

/* Filter visible citations  */
function filterCitations() {
   var rows = document.getElementById("citations").rows
   var year_min = parseInt($("year_min").value)
   var year_max = parseInt($("year_max").value)
   var title_filt = $("title_regex").value
   var abstract_filt = $("title_abstract_regex").value
   var exclude_filt = $("exclude_regex").value
   var journal_filt = $("journal_regex").value
   var r_title = RegExp(title_filt, "i")
   var r_abstract = RegExp(abstract_filt, "i")
   var r_exclude = RegExp(exclude_filt, "i")
   var r_journal = RegExp(journal_filt, "i")
   // row 0 is the heading
   for (var i = 1; i < rows.length; i++) {
      if (rows[i].className == "main" && rows[i].style.display != "none") {
         var cells = rows[i].cells;
         var title = cells[7].innerHTML;
         var abst = title;
         if (rows[i+2].length > 0) {
            abst += rows[i+2].cells[0].innerHTML; 
         }
         var year = parseInt(cells[4].innerHTML);
         var journal = cells[8].innerHTML;
         /* Apply filter criteria */
         if ( (year >= year_min && year <= year_max)
            && (title_filt == "" || r_title.test(title))
            && (abstract_filt == "" || r_abstract.test(abst))
            && (exclude_filt == "" || !r_exclude.test(abst))
            && (journal_filt == "" || r_journal.test(journal)))
         {
            // do nothing
         } else {
            rows[i].style.display = "none";            
         }
         rows[i+1].style.display = "none";
         rows[i+2].style.display = "none";            
      }
   }
}

/* Event handling for the citation filter */
function initCitationFilter() {
   $('clear_filter').onclick = hideTableRows;
   $('do_filter').onclick = filterCitations;
   $('filter_form').onsubmit = filterCitations;
}

/************************** WINDOW LOADING ******************************/

/* Hide the warning that scripts are blocked */
function hideScriptWarning() {
   $("script_warning").style.display = "none";
}

window.onload = function() {
   hideScriptWarning()
   make_toggles()
   initPubMedOpening()
   initTagSaving()
   initResultAppending()
   initTableEvents()
   initCitationFilter()
}
