import { IDsToResources } from '../../../utils';
import { ComponentWithResources, CompoundList } from '../../../index';
import ApiResourcePaginator from '../../ApiResourcePaginator';
import React from 'react';

export default function CompoundListFromAPI(props) {
  const molset = props.molset;

  const resourcesDef = IDsToResources(props.apiUrls.activitySetsRoot, molset.activities);
  const url = new URL(`${props.molset.id}/molecules/`, props.apiUrls.compoundSetsRoot);
  return (
    <ComponentWithResources
      {...props}
      definition={resourcesDef}
    >
      {
        (allLoaded, activitySets) => {
          return allLoaded ? (
            <ApiResourcePaginator
              url={url}
            >
              {
                currentItems => {
                  return (
                    <CompoundList
                      {...props}
                      activitySets={activitySets}
                      mols={currentItems}
                      paginate={false}
                    />
                  )
                }
              }
            </ApiResourcePaginator>
          ) : <div>Loading...</div>
        }
      }
    </ComponentWithResources>
  )
}