// graft() function
// Originally by Sean M. Burke from interglacial.com
// Closure support added by Maciek Adwent

function graft (parent, t, doc) {

    // Usage: graft( somenode, [ "I like ", ['em',
    //               { 'class':"stuff" },"stuff"], " oboy!"] )

    doc = (doc || parent.ownerDocument || document);
    var e;

    if(t == undefined) {
        throw complaining( "Can't graft an undefined value");
    } else if(t.constructor == String) {
        e = doc.createTextNode( t );
    } else if(t.length == 0) {
        e = doc.createElement( "span" );
        e.setAttribute( "class", "fromEmptyLOL" );
    } else {
        for(var i = 0; i < t.length; i++) {
            if( i == 0 && t[i].constructor == String ) {
                var snared;
                snared = t[i].match( /^([a-z][a-z0-9]*)\.([^\s\.]+)$/i );
                if( snared ) {
                    e = doc.createElement(   snared[1] );
                    e.setAttribute( 'class', snared[2] );
                    continue;
                }
                snared = t[i].match( /^([a-z][a-z0-9]*)$/i );
                if( snared ) {
                    e = doc.createElement( snared[1] );  // but no class
                    continue;
                }

                // Otherwise:
                e = doc.createElement( "span" );
                e.setAttribute( "class", "namelessFromLOL" );
            }

            if( t[i] == undefined ) {
                throw complaining("Can't graft an undefined value in a list! ");
            } else if(  t[i].constructor == String ||
                                    t[i].constructor == Array ) {
                graft( e, t[i], doc );
            } else if(  t[i].constructor == Number ) {
                graft( e, t[i].toString(), doc );
            } else if(  t[i].constructor == Object ) {
                // hash's properties => element's attributes
                for(var k in t[i]) {
                    // support for attaching closures to DOM objects
                    if(typeof(t[i][k])=='function'){
                        e[k] = t[i][k];
                    } else {
                        e.setAttribute( k, t[i][k] );
                    }
                }
            } else {
                throw complaining( "Object " + t[i] +
                    " is inscrutable as an graft arglet." );
            }
        }
    }

    parent.appendChild( e );
    return e; // return the topmost created node
}

function complaining (s) { alert(s); return new Error(s); }