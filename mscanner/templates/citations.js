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

/* Return whether event was a right click */
function is_rightclick(event) {
   var rightclick;
   if (event.which) {
      rightclick = (event.which == 3);
   }
   else if (event.button) {
      rightclick = (event.button == 2);
   }
   return rightclick;
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

// 0 = Unclassified (gray)
// 1 = Negative (blue)
// 2 = Positive (red)

// Read classifications from disk
function loadTags(fname) 
{
   var text = readFile(fname)
   if (text == null) {
      alert("Failed to read text from file");
      return;
   }
   var lines = text.split("\n");
   for (i = 0; i < lines.length; i++) {
      var values = lines[i].split(",");
      var pmid = values[0];
      var score = values[1];
      var classification = values[2];
      switch(classification) 
      {
      case " ": setTag(pmid, 0); break;
      case "0": setTag(pmid, 1); break;
      case "1": setTag(pmid, 2); break;
      }
   }
   alert("Successfully loaded classifications");
}

// Save classifications to disk
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
   fname = document.getElementById("fname").value;
   if (writeFile(fname, text) == true) {
      alert("Successfully saved classifications");
   }
}

/* 
Read citations from disk and append to the table
If append is false, we prepend the data to the table instead.
*/
function loadCitations(fname, append) {
   var text = readFile(fname);
   var table = text.match(/<[t]body>[^\v]*?<\/tbody>/);
   var blob = table[0].substr(7,table[0].length-15);
   var tbody = document.getElementById("citations").getElementsByTagName("tbody")[0]
   alert(tbody.innerHTML);
   if (append == true) {
      tbody.innerHTML = tbody.innerHTML + blob;
   } else {
      tbody.innerHTML = blob + tbody.innerHTML;
   }
   /* because new text is unhidden, and no events are linked */
   hide_allrows();
   add_table_events();
}

/* Open positive documents in PubMed */
function openPubMed(all) {
   var rows = document.getElementById("citations").rows;
   var qstring = "http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&amp;db=pubmed&amp;list_uids=";
   var count = 0;
   for (var i = 0; i < rows.length; i++) {
      if (rows[i].className == "main") {
         var pmid = rows[i].id;
         if (getTag(pmid) == 2 || all == true) {
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

// Get classification tag (0,1,2)
function getTag(pmid) {
   var color = document.getElementById(pmid).cells[0].style.backgroundColor
   switch (color) {
      case "": return 0;
      case "blue": return 1;
      case "red": return 2;
      default: throw "Invalid tag color " + color;
   } 
}

// Set classification tag (takes 0,1,2)
function setTag(pmid, value) {
   var row = document.getElementById(pmid)
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

/* Hides rows not matching the given regex filter. */
function filterCitations(filter) {
   var rows = document.getElementById("citations").rows;
   // row 0 is the heading
   for (var i = 1; i < rows.length; i++) {
      if (rows[i].className == "main") {
         var cells = rows[i].cells;
         var text = cells[7].innerHTML
         if (rows[i+2].length > 0) {
            text += rows[i+2].cells[0].innerHTML;
         }
         if (filter == "" || RegExp(filter, "i").test(text)) {
            rows[i].style.display = "";            
         } else {
            rows[i].style.display = "none";            
         }
         rows[i+1].style.display = "none";
         rows[i+2].style.display = "none";            
      }
   }
}

/* Hides rows that are outside the given range of years */
function filterByYear(startyear, endyear) {
   var rows = document.getElementById("citations").rows;
   var start = parseInt(startyear);
   var finish = parseInt(endyear);
   // row 0 is the heading
   for (var i = 1; i < rows.length; i++) {
      if (rows[i].className == "main") {
         var cells = rows[i].cells;
         var year = parseInt(cells[4].innerHTML);
         if (year >= startyear && year <= endyear) {
            rows[i].style.display = "";            
         } else {
            rows[i].style.display = "none";            
         }
         rows[i+1].style.display = "none";
         rows[i+2].style.display = "none";            
      }
   }
   return false;
}

/* Hide all expanded author/abstract rows in the citations table */
function hide_allrows() {
   var rows = document.getElementById("citations").rows;
   // row 0 is the heading
   for (i = 1; i < rows.length; i+=3) {
      rows[i+1].style.display = "none";            
      rows[i+2].style.display = "none";
   }
}

/* Adds the events for classification/author/abstract display */
function add_table_events() {

   /* Event handler for classification cells.
   Cycles forward on left clicks, backward on right clicks */
   function onclick_classify(e) {
      if (!e) var e = window.event;
      // Assuming target is cell, with parent row, with PMID as id
      var pmid = this.parentNode.id;
      if (is_rightclick(e)) {
         setTag(pmid, (getTag(pmid)+1) % 3);
      } else {
         setTag(pmid, (getTag(pmid)+2) % 3);
      }
      return false; // suppress default click action
   }

   /* Event handler for author-toggle */
   function onclick_author() {
      // Author row is sibling of target's parent row.
      toggle(this.parentNode.nextSibling);
   }

   /* Event handler for abstract-toggle */
   function onclick_abstract() {
      // Abstract row is second sibling of target's parent row.
      toggle(this.parentNode.nextSibling.nextSibling);
   }
   
   var rows = document.getElementById("citations").rows;
   // row 0 is the heading
   for (var i = 1; i < rows.length; i+=3) {
      rows[i].cells[0].onclick = onclick_classify;
      rows[i].cells[5].onclick = onclick_author;
      rows[i].cells[6].onclick = onclick_abstract;
   }
}

/* add event handlers to citation form */
function add_form_events() {
   $('save_tags').onclick = function() {
      saveTags($('fname').value);
   }
   $('load_tags').onclick = function() {
      loadTags($('fname').value);
   }
   $('pubmed_relevant').onclick = function() {
      openPubMed(false);
   }
   $('pubmed_all').onclick = function() {
      openPubMed(true);
   }
   $('append_html').onclick = function() {
      loadCitations($('fname').value);
   }
   $('clear_filter').onclick = function() {
      $('filter').value = '';
      filterCitations('');
   }
   $('do_filter').onclick = function() {
      filterCitations($('filter').value);
   }
   $('do_year_filter').onclick = function() {
      filterByYear($('yearstart').value, $('yearend').value);
   }
}

/* things to do once the page is loaded */
window.onload = function() {
   add_table_events();
   add_form_events();
   hide_allrows();
}
