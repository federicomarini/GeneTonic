shinyjs.loadStringData = function(params) {

    var defaultParams = {
        organism : '9606',
        gene : 'ACTB'
    };

    params = shinyjs.getParams(params, defaultParams);

    getSTRING('https://string-db.org', {
        'species': params.organism,
        'identifiers': [params.gene],
        'network_flavor':'confidence'
    });
};
