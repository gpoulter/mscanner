function toggle(id, type) {
   if (document.getElementById(id).style.display == "none") {
      show(id, type);
   } else {
      hide(id);
   }
}

function hide(id) {
   document.getElementById(id).style.display = "none";
}

function show(id, type) {
   if (type == null) {
      document.getElementById(id).style.display = "";
   } else {
      document.getElementById(id).style.display = type;
   }
}  
