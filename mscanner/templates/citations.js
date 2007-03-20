// 0 = Unclassified (gray)
// 1 = Negative (blue)
// 2 = Positive (red)

// Read classifications from disk
function loadTags() {
   fname = document.getElementById("fname").value
   text = readFile(fname);
   lines = text.split("\n");
   for(i = 0; i < lines.length; i++) {
      values = lines[i].split(",");
      pmid = values[0];
      score = values[1];
      classification = values[2];
      switch(classification) {
      case " ": setTag(pmid, 0); break;
      case "0": setTag(pmid, 1); break;
      case "1": setTag(pmid, 2); break;
      }
   }
}

// Save classifications to disk
function saveTags() {
   text = "";
   rows = document.getElementById("citations").rows;
   for (i = 0; i < rows.length; i++) {
      if (rows[i].className == "main") {
         pmid = rows[i].id;
         score = document.getElementById("s_"+pmid).innerHTML;
         line = pmid+","+score+","
         switch(getTag(pmid)) {
         case 0: line += " \n"; break;
         case 1: line += "0\n"; break;
         case 2: line += "1\n"; break;
         }
         text += line;
      }
   }
   fname = document.getElementById("fname").value
   writeFile(fname, text);
}

// Open positive documents in PubMed
function openPubMed() {
   rows = document.getElementById("citations").rows;
   qstring = "http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&amp;db=pubmed&amp;list_uids=";
   count = 0;
   for (i = 0; i < rows.length; i++) {
      if (rows[i].className == "main") {
         pmid = rows[i].id;
         if (getTag(pmid) == 2) {
            if (count > 0) {
               qstring += ",";
            }
            qstring += pmid;
            count++;
         }
      }
   }
   window.open(qstring);
}

// Get classification tag (0,1,2)
function getTag(pmid) {
   color = document.getElementById("c_"+pmid).style.backgroundColor
   switch(color) {
   case "": return 0;
   case "blue": return 1;
   case "red": return 2;
   default: throw "Invalid tag color " + color;
   } 
}

// Set classification tag (takes 0,1,2)
function setTag(pmid, value) {
   tagstyle = document.getElementById("c_"+pmid).style
   rowstyle = document.getElementById(pmid).style
   switch(value) {
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

// Cycle the tag through unclassified/negative/positive
function classify(pmid) {
   setTag(pmid, (getTag(pmid)+1) % 3);
}

// Filter citations (empty string shows all citations)
function filterCitations(filter) {
   rows = document.getElementById("citations").rows;
   for (i = 0; i < rows.length; i++) {
      if (rows[i].className == "main") {
         pmid = rows[i].id;
         if (filter == "" || 
            RegExp(filter, "i").test($("t_"+pmid).innerHTML + $("ab_"+pmid).innerHTML)
         ) {
            rows[i].style.display = "";            
         } else {
            rows[i].style.display = "none";            
         }
         rows[i+1].style.display = "none";            
         rows[i+2].style.display = "none";            
      }
   }
}

