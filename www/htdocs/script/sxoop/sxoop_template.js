var SXOOP = SXOOP || {};

SXOOP.template = {};

/**
 * Does this need explaining  ?
 */
SXOOP.template.map = function(array, func)
{
   var result = [];
   for (var i = 0;i < array.length; i++){ 
       var mapped = func(array[i]);
       for (var j = 0; j < mapped.length; j++){
           result.push(mapped[j]);
       }
   }
   return result;
};

/**
 * Parse a template (supplied as a string), substituting
 * the supplied object ($_) 
 * The $_ variable refers to the object which was passed into the parse function
 * Of course, all other global variables/functions are accessible too.
 */
SXOOP.template.parse = function(str,$_)
{
    var singleLine = str.replace(/[\n\r]/g,"");
    // innerHTML automatically converts < to &lt; and > to &gt;
    singleLine = singleLine.replace(/&lt;/g,"<");
    singleLine = singleLine.replace(/&gt;/g,">");

    /**
     * The include function facilitates inclusion of inner templates
     * Note: This include function is local to the parse function and will
     * override a global include function in the template scope only.
     */
    var include = function(elementId){
        var included = document.getElementById(elementId).innerHTML;
        return SXOOP.template.parse(included,$_);
    };
   
    /**
     * Split the template into parts
     */
    var parts = SXOOP.template.map(singleLine.split("[:"),function (part){
        var result = [];
        if (part.match(/:\]/)){
            result = part.split(/:\]/g);
            result[0] = "[:" + result[0];
        }else{
            result = [part];
        }
        return result;
    });
    /**
     * In firefox the following would suffice instead.
     * IE's implementation of split() is broken -  doesn't retain captured parts.
     *
     * parts = singleLine.split(/(\[:.*?):\]/);
     *
     * Process each part
     */
    var result = SXOOP.template.map(parts,function (part){
        var result = "";
        if (part.match(/\[:=/)){
            var inner = part.replace(/^\[:=\s*/,"");
            return ["theArray.push(" + inner + ");"];
        }
        if (part.match(/^\[:/)){
            var inner = part.replace(/^\[:/,"");
            return [inner];
        }else{
            part = part.replace(/\"/g,"\\\"");
            return ["theArray.push(\"" + part + "\");"];
        }
    });
    var theArray = [];
    result.push("theArray.join('');");
    var javascript = result.join("\n");
    return eval(javascript);
};

