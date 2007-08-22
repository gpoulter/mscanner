/* helper shortcut for getElementById */
function $(name) {
   return document.getElementById(name);
}

/* Return target of the event */
function event_target(event) {
   var targ;
   if (event.target) {
      targ = event.target; // Netscape
   }
   else if (event.srcElement) {
      targ = event.srcElement; // Microsoft
   }
   if (targ.nodeType == 3) { // defeat Safari bug
      targ = targ.parentNode;
   }
   return targ;
}  

visible_element = null;

/* toggle element */
function toggle(element, vis_state) {
   if (element.style.display == "none") {
      show(element, vis_state);
   } else {
      hide(element);
   }
}

/* hide element, setting last_visible to null */
function hide(element) {
   element.style.display = "none";
   if (visible_element == element) {
      visible_element = null;
   }
}

/* unhide specified element, and assign it to last_visible.
If last_visible is non-null, we hide it */
function show(element, vis_state) {
   if (visible_element != null && visible_element != element) {
      hide(visible_element);
   }
   if (vis_state == null) {
      element.style.display = "";
   } else {
      element.style.display = vis_state;
   }
   visible_element = element;
}

/* given an event, return its target in a cross-browser manner */
function getEventTarget(e) {
   var targ;
   if (!e) var e = window.event;
   if (e.target) targ = e.target;
   else if (e.srcElement) targ = e.srcElement;
   if (targ.nodeType == 3) // defeat Safari bug
      targ = targ.parentNode;
   return targ;
}

/* create events for all toggle buttons */
function make_toggles() {
  /* handler for toggle buttons, using the id's like toggle_x */
  function toggle_handler() {
    toggle($(this.id.substr(7)));
  }
  /* find the toggle buttons and link them up*/
  var buttons = document.getElementsByTagName("button");
  for(var i = 0; i < buttons.length; i++) {
    button = buttons[i];
    if(button.className == "toggle") {
      target = $(button.id.substr(7));
      if(target) {
         button.onclick = toggle_handler;
         /* initially hide the thing being toggled */
         hide(target);
      }
    }
  }
}
