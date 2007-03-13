function toggle(id) {
   if (document.getElementById(id).style.display == "none") {
      show(id);
   } else {
      hide(id);
   }
}
function hide(id) {
   document.getElementById(id).style.display = "none";
}
function show(id) {
   document.getElementById(id).style.display = "";
}  
