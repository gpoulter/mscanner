// Return element for ID
function $(name) {
   return document.getElementById(name);
}

visible_id = null

// Toggle such that only one element stays visible
function toggle(id, vis_state) {
   if ($(id).style.display == "none") {
      show(id, vis_state);
   } else {
      hide(id);
   }
}

// Hide element, clearing last_visible
function hide(id) {
   $(id).style.display = "none";
   if (visible_id == id) {
      visible_id = null;
   }
}

// Show element, changing last_visible
function show(id, vis_state) {
   if (visible_id != null && visible_id != id) {
      hide(visible_id);
   }
   if (vis_state == null) {
      $(id).style.display = "";
   } else {
      $(id).style.display = vis_state;
   }
   visible_id = id;
}
