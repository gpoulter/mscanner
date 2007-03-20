function $(name) {
   return document.getElementById(name);
}

function toggle(id, type) {
   if ($(id).style.display == "none") {
      show(id, type);
   } else {
      hide(id);
   }
}

function hide(id) {
   $(id).style.display = "none";
}

function show(id, type) {
   if (type == null) {
      $(id).style.display = "";
   } else {
      $(id).style.display = type;
   }
}

